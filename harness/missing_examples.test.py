"""Controls and fixed-seed calibration for category G (missing data and bias).

multiple-imputation-chained-equations (R only), delta-adjusted-tipping-point, e-value-quantitative-bias,
inverse-probability-censoring-weights. Each fixture regenerates byte for byte and each truth is
recomputed. Oracles: Rubin's rules and the Barnard-Rubin degrees of freedom recomputed from mice's
pooled components; the tipping point sits exactly where the lower limit is zero, and the effect moves
linearly in delta with slope equal to the treated arm's missing fraction; the E-value formula reproduces
the published worked example (VanderWeele and Ding 2017: RR 3.9, lower limit 1.8, E-values 7.26 and 3.0);
the IPCW sandwich is rebuilt with numerical derivatives. Negative controls over the seeds: each entry's
naive analysis is biased.

SEEDS: 50 by default, as CI runs it; MALLARD_SEEDS=100 for the full calibration each expected.json
records (owner decision E5, 2026-09-23). Run with --measure to print the numbers. Locally, MALLARD_R_LIBS
may name an R library to put first on the path.
"""
import contextlib
import csv
import importlib.util
import io
import json
import math
import os
import runpy
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy import integrate, stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check import parse_harness_block  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
SEEDS = range(int(os.environ.get("MALLARD_SEEDS", "50")))


def rank_quantile(misses):
    """The recovery statistic: the quantile at the same ORDER STATISTIC as the 99th percentile of 100 draws.

    np.quantile(q=0.99) over 100 misses lands just above the second-largest; over 50 it lands halfway to
    the largest, which is stricter than the calibration the tolerances were set from. The level
    (n - 1.99) / (n - 1) is exactly 0.99 at n = 100 and the same position at any n (owner decision E6,
    2026-09-23).
    """
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


def run_r_many(name, frames):
    libs = os.environ.get("MALLARD_R_LIBS")
    prelude = f'.libPaths(c("{libs}", .libPaths())); ' if libs else ""
    with tempfile.TemporaryDirectory() as tmp:
        for i, generated in enumerate(frames):
            Path(tmp, str(i)).mkdir()
            Path(tmp, str(i), "fixture.csv").write_text(csv_text(generated))
        script = (ROOT / "lib" / name / "r.R").as_posix()
        code = (prelude + f'for (i in 0:{len(frames) - 1}) {{ setwd(file.path("{tmp}", i)); '
                f'out <- capture.output(source("{script}")); cat("@@", i, "\\n"); cat(out, sep = "\\n") }}')
        text = subprocess.check_output(["Rscript", "-e", code], text=True, stderr=subprocess.DEVNULL)
    return [parse_harness_block(chunk) for chunk in text.split("@@")[1:]]


def miss99(values, truth):
    return rank_quantile([abs(v - truth) for v in values])


def normal_mean(f):
    return integrate.quad(lambda x: f(x) * math.exp(-0.5 * x * x) / math.sqrt(2 * math.pi), -12, 12,
                          epsabs=1e-14, epsrel=1e-13, limit=200)[0]


# ---------------------------------------------------------------------------------------------------

def multiple_imputation(measured, truths):
    name = "multiple-imputation-chained-equations"
    regenerates(name)
    m = load(name)
    truths[name] = {"effect": m.EFFECT}
    runs = run_r_many(name, [list(m.rows())] + [list(m.rows(s)) for s in SEEDS])
    out, runs = runs[0], runs[1:]
    # Rubin's rules and Barnard-Rubin df, recomputed from the pooled components.
    mm, ubar, b = out["imputations"], out["within_variance"], out["between_variance"]
    t = ubar + (1 + 1 / mm) * b
    assert abs(math.sqrt(t) - out["effect_se"]) < 1e-9
    lam = (1 + 1 / mm) * b / t
    nu_com = out["n"] - 3
    nu_old = (mm - 1) / lam ** 2
    nu_obs = (nu_com + 1) / (nu_com + 3) * nu_com * (1 - lam)
    df = nu_old * nu_obs / (nu_old + nu_obs)
    assert abs(df / out["df"] - 1) < 1e-7, f"Barnard-Rubin df {df} != mice {out['df']}"  # relative: ubar and b are printed to 10 dp
    q = stats.t.ppf(0.975, df)
    assert abs(out["effect"] - q * math.sqrt(t) - out["effect_lcl"]) < 1e-8
    assert np.mean([r["complete_case_effect"] for r in runs]) - m.EFFECT < -0.1, "complete-case should be biased"
    measured[name] = {"effect": miss99([r["effect"] for r in runs], m.EFFECT),
                      "complete_case_mean": float(np.mean([r["complete_case_effect"] for r in runs])),
                      "mi_mean": float(np.mean([r["effect"] for r in runs]))}


def tipping_point(measured, truths):
    name = "delta-adjusted-tipping-point"
    regenerates(name)
    m = load(name)
    truths[name] = {"effect": m.EFFECT}
    out, scope = run_python(name)
    effect_at, z = scope["effect_at"], scope["z"]
    tip = out["tipping_delta"]
    e, s = effect_at(tip)
    assert abs(e - z * s) < 1e-8, "the lower limit is not zero at the tipping point"
    assert effect_at(tip + 0.1)[0] - z * effect_at(tip + 0.1)[1] > 0 > effect_at(tip - 0.1)[0] - z * effect_at(tip - 0.1)[1]
    d = scope["d"]
    frac = float((d.loc[d.treated == 1, "r"] == 0).mean())
    assert abs((effect_at(-1.0)[0] - effect_at(0.0)[0]) - (-frac)) < 1e-10, "effect should move by the missing fraction"
    runs = [run_python(name, list(m.rows(s)))[0] for s in SEEDS]
    assert np.mean([r["complete_case_difference"] for r in runs]) - m.EFFECT < -0.2, "complete-case should be biased"
    measured[name] = {"effect": miss99([r["effect"] for r in runs], m.EFFECT)}


def e_value(measured, truths):
    name = "e-value-quantitative-bias"
    regenerates(name)
    m = load(name)
    truths[name] = {"log_rr_bias_adjusted": math.log(m.RR_A)}
    out, scope = run_python(name)
    ev = scope["e_value"]
    assert round(ev(3.9), 2) == 7.26 and round(ev(1.8), 1) == 3.0, "published worked example"
    assert abs(ev(0.8) - (1.25 + math.sqrt(1.25 * 0.25))) < 1e-12, "a protective RR is inverted"
    b = m.bias_parameters()
    assert abs(b["p1"] - scope["p1"]) < 5e-5 and abs(b["p0"] - scope["p0"]) < 5e-5 and b["rr_uy"] == scope["rr_uy"]
    runs = [run_python(name, list(m.rows(s)))[0] for s in SEEDS]
    assert np.mean([r["log_rr_observed"] for r in runs]) - math.log(m.RR_A) > 0.15, "the crude RR should be confounded"
    measured[name] = {"log_rr_bias_adjusted": miss99([r["log_rr_bias_adjusted"] for r in runs], math.log(m.RR_A))}


def ipcw(measured, truths):
    name = "inverse-probability-censoring-weights"
    regenerates(name)
    m = load(name)
    rd = normal_mean(lambda l: m.risk(1, l)) - normal_mean(lambda l: m.risk(0, l))
    assert abs(m.truth()["risk_difference"] - rd) < 1e-10
    truths[name] = {"risk_difference": rd}
    out, scope = run_python(name)
    X, a, r, y = scope["X"], scope["a"], scope["r"], scope["y"]
    theta = np.r_[np.asarray(scope["cm"].params), scope["mu1"], scope["mu0"]]
    k = X.shape[1]

    def psi(th):
        p = 1 / (1 + np.exp(-X @ th[:k]))
        w = r / p
        return np.column_stack([X * (r - p)[:, None], w * a * (y - th[k]), w * (1 - a) * (y - th[k + 1])])

    A = np.column_stack([(psi(theta + h).sum(0) - psi(theta - h).sum(0)) / 2e-6 for h in np.eye(len(theta)) * 1e-6])
    Ai = np.linalg.inv(A)
    V = Ai @ (psi(theta).T @ psi(theta)) @ Ai.T
    c = np.r_[np.zeros(k), 1.0, -1.0]
    assert abs(math.sqrt(c @ V @ c) - out["risk_difference_se"]) < 1e-7
    runs = [run_python(name, list(m.rows(s)))[0] for s in SEEDS]
    assert np.mean([x["complete_case_risk_difference"] for x in runs]) - rd > 0.03, "complete-case should be biased"
    measured[name] = {"risk_difference": miss99([x["risk_difference"] for x in runs], rd)}


def check_tolerances(measured, truths):
    for name, truth in truths.items():
        exp = json.loads((ROOT / "lib" / name / "expected.json").read_text())
        for key, value in truth.items():
            assert abs(exp["truth"][key] - value) < 1e-10, f"{name} {key}: expected.json {exp['truth'][key]} != {value}"
            tol = float(exp["recovery"].get("tolerance_by_key", {}).get(key, exp["recovery"]["tolerance_estimate"]))
            assert measured[name][key] < tol, f"{name} {key}: 99th percentile miss {measured[name][key]:.4f} >= {tol}"


if __name__ == "__main__":
    measured, truths = {}, {}
    multiple_imputation(measured, truths)
    tipping_point(measured, truths)
    e_value(measured, truths)
    ipcw(measured, truths)
    if "--measure" in sys.argv:
        print(json.dumps({"seeds": len(SEEDS), "measured": measured, "truths": truths}, indent=2))
        sys.exit(0)
    check_tolerances(measured, truths)
    print(json.dumps({"entries": len(measured), "seeds": len(SEEDS), "result": "pass"}, indent=2))
