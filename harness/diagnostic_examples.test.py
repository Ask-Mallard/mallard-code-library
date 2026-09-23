"""Controls and fixed-seed calibration for batch H-1 (diagnostic accuracy and calibration).

roc-auc-delong, comparing-two-diagnostic-tests, predictive-values-likelihood-ratios,
calibration-slope-intercept. Each fixture regenerates byte for byte and each truth is recomputed
(calibration by scipy quadrature). Oracles: the AUC equals scipy's Mann-Whitney U / (m n), and the
reversed direction gives 1 - AUC; McNemar's statistic equals (b - c)^2 / (b + c); the Wilson interval and
the likelihood-ratio SEs are rebuilt from their formulas; the PPV at another prevalence follows Bayes'
rule from the estimated sensitivity and specificity. Negative controls over the seeds: the Youden index
at the threshold chosen on the data is optimistic; the paired SE of the AUC difference is below the
unpaired one; the intercept of the slope model is not calibration-in-the-large.

SEEDS: 50 by default, as CI runs it; MALLARD_SEEDS=100 for the recorded calibration (owner decisions E5,
E6). Run with --measure to print the numbers.
"""
import contextlib
import csv
import importlib.util
import io
import json
import math
import os
import runpy
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy import integrate, optimize, stats
from scipy.special import expit

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check import parse_harness_block  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
SEEDS = range(int(os.environ.get("MALLARD_SEEDS", "50")))  # 50 in CI; MALLARD_SEEDS=100 reproduces the recorded calibration (owner decisions E5, E6)


def rank_quantile(misses):
    """The quantile at the same order statistic as the 99th percentile of 100 draws (owner decision E6)."""
    n = len(misses)
    return float(np.quantile(misses, (n - 1.99) / (n - 1)))


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "lib" / name / "fixture.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def csv_text(generated):
    out = io.StringIO(newline="")
    writer = csv.DictWriter(out, fieldnames=list(generated[0]))
    writer.writeheader()
    writer.writerows(generated)
    return out.getvalue()


def regenerates(name):
    assert csv_text(list(load(name).rows())).encode() == (ROOT / "lib" / name / "fixture.csv").read_bytes(), \
        f"{name}: fixture.csv does not regenerate byte for byte"


def run_python(name, generated=None):
    script = ROOT / "lib" / name / "python.py"
    with tempfile.TemporaryDirectory() as tmp:
        src = (ROOT / "lib" / name / "fixture.csv").read_text() if generated is None else csv_text(generated)
        Path(tmp, "fixture.csv").write_text(src)
        here = os.getcwd()
        buf = io.StringIO()
        try:
            os.chdir(tmp)
            with contextlib.redirect_stdout(buf):
                scope = runpy.run_path(str(script))
        finally:
            os.chdir(here)
    return parse_harness_block(buf.getvalue()), scope


def runs(name):
    m = load(name)
    return [run_python(name, list(m.rows(s)))[0] for s in SEEDS]


def recovery(results, truth):
    return {key: rank_quantile([abs(r[key] - value) for r in results]) for key, value in truth.items()}


# ---------------------------------------------------------------------------------------------------

def roc(measured, truths):
    name = "roc-auc-delong"
    regenerates(name)
    m = load(name)
    auc = stats.norm.cdf(m.SHIFT / math.sqrt(2))
    assert abs(auc - m.true_auc()) < 1e-12
    truths[name] = {"auc": auc}
    out, scope = run_python(name)
    x, y = scope["x"], scope["y"]
    u = stats.mannwhitneyu(x, y, alternative="two-sided").statistic
    assert abs(u / (len(x) * len(y)) - out["auc"]) < 1e-9, "AUC is not the Mann-Whitney probability"
    reversed_psi = (y[None, :] > x[:, None]) + 0.5 * (x[:, None] == y[None, :])
    assert abs(reversed_psi.mean() - (1 - out["auc"])) < 1e-9
    r = runs(name)
    best_youden = 2 * stats.norm.cdf(m.SHIFT / 2) - 1
    optimism = np.mean([x["sensitivity"] + x["specificity"] - 1 for x in r]) - best_youden
    assert optimism > 0.01, "the Youden index at a threshold chosen on the data should be optimistic"
    measured[name] = recovery(r, truths[name])
    measured[name]["youden_optimism"] = float(optimism)


def two_tests(measured, truths):
    name = "comparing-two-diagnostic-tests"
    regenerates(name)
    m = load(name)
    t = m.truth()
    truths[name] = {"auc_difference": t["auc_difference"], "sensitivity_difference": t["sensitivity_difference"]}
    assert abs(t["auc_difference"] - (stats.norm.cdf(1.2 / math.sqrt(2)) - stats.norm.cdf(0.7 / math.sqrt(2)))) < 1e-12
    out, _ = run_python(name)
    b, c = out["discordant_a_only"], out["discordant_b_only"]
    assert abs(out["mcnemar_statistic"] - (b - c) ** 2 / (b + c)) < 1e-10
    r = runs(name)
    assert all(x["auc_difference_se"] < x["auc_difference_se_unpaired"] for x in r), "pairing should shrink the SE"
    measured[name] = recovery(r, truths[name])


def predictive(measured, truths):
    name = "predictive-values-likelihood-ratios"
    regenerates(name)
    m = load(name)
    truths[name] = m.truth()
    assert abs(truths[name]["ppv"] - 0.68) < 1e-12 and abs(truths[name]["npv"] - 0.96) < 1e-12
    out, s = run_python(name)
    z = stats.norm.ppf(0.975)

    def wilson(x, n):
        p = x / n
        centre = (p + z * z / (2 * n)) / (1 + z * z / n)
        half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
        return centre - half, centre + half

    lo, hi = wilson(s["tp"], s["tp"] + s["fp"])
    assert abs(lo - out["ppv_lcl"]) < 1e-10 and abs(hi - out["ppv_ucl"]) < 1e-10
    se, sp = s["tp"] / s["nd"], s["tn"] / s["nnd"]
    assert abs(math.sqrt((1 - se) / (se * s["nd"]) + sp / ((1 - sp) * s["nnd"])) - out["log_lr_positive_se"]) < 1e-9
    other = m.truth(0.05)["ppv"]
    assert abs(other - 0.85 * 0.05 / (0.85 * 0.05 + 0.10 * 0.95)) < 1e-12
    assert abs(out["ppv_at_prevalence_0.05"] - se * 0.05 / (se * 0.05 + (1 - sp) * 0.95)) < 1e-9
    measured[name] = recovery(runs(name), truths[name])


def calibration(measured, truths):
    name = "calibration-slope-intercept"
    regenerates(name)
    m = load(name)
    b = m.TRUE[1] / m.PUBLISHED[1]
    a = m.TRUE[0] - b * m.PUBLISHED[0]
    sd = math.hypot(m.PUBLISHED[1], m.PUBLISHED[2])
    mean = lambda f: integrate.quad(lambda v: f(v) * stats.norm.pdf(v, m.PUBLISHED[0], sd),
                                    m.PUBLISHED[0] - 12 * sd, m.PUBLISHED[0] + 12 * sd, epsabs=1e-14, limit=200)[0]
    observed = mean(lambda v: expit(a + b * v))
    citl = optimize.brentq(lambda c: mean(lambda v: expit(c + v)) - observed, -5, 5, xtol=1e-14)
    truth = {"calibration_slope": b, "calibration_in_the_large": citl, "observed_expected": observed / mean(expit)}
    stated = m.truth()
    for key, value in truth.items():
        assert abs(stated[key] - value) < 1e-9, (key, stated[key], value)
    truths[name] = truth
    r = runs(name)
    assert np.mean([abs(x["slope_model_intercept"] - x["calibration_in_the_large"]) for x in r]) > 0.05, \
        "the slope model's intercept should not be calibration-in-the-large"
    measured[name] = recovery(r, truth)


def check_tolerances(measured, truths):
    for name, truth in truths.items():
        exp = json.loads((ROOT / "lib" / name / "expected.json").read_text())
        for key, value in truth.items():
            assert abs(exp["truth"][key] - value) < 1e-9, f"{name} {key}: expected.json {exp['truth'][key]} != {value}"
            tol = float(exp["recovery"].get("tolerance_by_key", {}).get(key, exp["recovery"]["tolerance_estimate"]))
            assert measured[name][key] < tol, f"{name} {key}: recovery miss {measured[name][key]:.4f} >= {tol}"


if __name__ == "__main__":
    measured, truths = {}, {}
    roc(measured, truths)
    two_tests(measured, truths)
    predictive(measured, truths)
    calibration(measured, truths)
    if "--measure" in sys.argv:
        print(json.dumps({"seeds": len(SEEDS), "measured": measured, "truths": truths}, indent=2))
        sys.exit(0)
    check_tolerances(measured, truths)
    print(json.dumps({"entries": len(measured), "seeds": len(SEEDS), "result": "pass"}, indent=2))
