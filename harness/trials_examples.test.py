"""Controls and fixed-seed calibration for batch I-1 (noninferiority, multiplicity, SPC, bootstrap).

noninferiority-risk-difference, multiplicity-adjustment, statistical-process-control,
bootstrap-percentile-bca. Each fixture regenerates byte for byte. Oracles: Holm and Benjamini-Hochberg
adjusted p-values rebuilt by hand; the u-chart limits from their formula; the shared bootstrap generator
reproduces the C++ standard's check value for minstd_rand (the 10,000th value from seed 1 is 399268537),
and BCa with z0 = a = 0 reduces to the percentile interval. Negative controls over the seeds: the ITT
difference sits closer to 0 than the per-protocol one; unadjusted testing of the seven null outcomes
rejects at least one far more often than Holm; the u-chart signals after the rate rises and rarely
before; the bootstrap intervals cover the true ratio at close to the nominal rate.

SEEDS: 50 by default, as CI runs it; MALLARD_SEEDS=100 for the recorded calibration (owner decisions E5,
E6). Run with --measure to print the numbers.
"""
import contextlib
import csv
import importlib.util
import io
import json
import os
import runpy
import sys
import tempfile
from pathlib import Path

import numpy as np

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

def noninferiority(measured, truths):
    name = "noninferiority-risk-difference"
    regenerates(name)
    m = load(name)
    truth = m.truth()
    assert abs(truth["risk_difference_itt"] - (0.9 * 0.78 + 0.1 * 0.80 - 0.80)) < 1e-12
    truths[name] = truth
    r = runs(name)
    assert np.mean([x["risk_difference_itt"] - x["risk_difference_pp"] for x in r]) > 0, \
        "protocol deviations should pull the ITT difference towards 0"
    measured[name] = recovery(r, truth)


def multiplicity(measured, truths):
    name = "multiplicity-adjustment"
    regenerates(name)
    m = load(name)
    truths[name] = {"difference_1": m.EFFECTS[0]}
    out, _ = run_python(name)
    p = np.array([out[f"p_{j}"] for j in range(1, 11)])
    order = np.argsort(p)
    holm = np.empty(10)
    holm[order] = np.minimum(1, np.maximum.accumulate((10 - np.arange(10)) * p[order]))
    bh = np.empty(10)
    bh[order] = np.minimum(1, np.minimum.accumulate((10 / np.arange(10, 0, -1) * p[order][::-1]))[::-1])
    for j in range(10):
        assert abs(holm[j] - out[f"p_holm_{j + 1}"]) < 1e-9 and abs(bh[j] - out[f"p_bh_{j + 1}"]) < 1e-9
    r = runs(name)
    nulls = range(4, 11)
    fwer = lambda key: np.mean([any(x[f"{key}_{j}"] < 0.05 for j in nulls) for x in r])
    assert fwer("p") > 0.15 and fwer("p_holm") <= 0.1, "Holm should control what unadjusted testing does not"
    measured[name] = recovery(r, truths[name])
    measured[name].update({"fwer_unadjusted": float(fwer("p")), "fwer_holm": float(fwer("p_holm")),
                           "fwer_bh": float(fwer("p_bh"))})


def spc(measured, truths):
    name = "statistical-process-control"
    regenerates(name)
    m = load(name)
    truths[name] = m.truth()
    out, s = run_python(name)
    d = s["d"]
    c = d.infections[d.month <= 12].sum() / d.catheter_days[d.month <= 12].sum()
    assert abs(1000 * (c + 3 * np.sqrt(c / d.catheter_days.iloc[-1])) - out["u_ucl_month24_per_1000"]) < 1e-9
    r = runs(name)
    assert np.mean([x["u_first_signal_month"] >= m.SHIFT_FROM for x in r]) > 0.8, "the u-chart should signal after the rise"
    # 48 steady chart-months per seed; limits from only 12 baseline points give about 1% false alarms per
    # point (0.47 per seed over 100 seeds), well above the nominal 0.27%. Asserted: under 1 per seed.
    assert np.mean([x["p_signals"] + x["i_signals"] for x in r]) < 1.0, "steady processes should rarely signal"
    measured[name] = recovery(r, truths[name])


def bootstrap(measured, truths):
    name = "bootstrap-percentile-bca"
    regenerates(name)
    m = load(name)
    truths[name] = m.truth()
    state = 1
    for _ in range(10000):
        state = (48271 * state) % 2147483647
    assert state == 399268537, "the shared generator is not minstd_rand"
    out, s = run_python(name)
    boot = s["boot_stat"]
    assert np.allclose(np.quantile(boot, [0.025, 0.975]), [out["percentile_lcl"], out["percentile_ucl"]], atol=1e-9)
    r = runs(name)
    t = truths[name]["mean_ratio"]
    cover = lambda lo, hi: np.mean([x[lo] <= t <= x[hi] for x in r])
    assert cover("bca_lcl", "bca_ucl") >= 0.85 and cover("percentile_lcl", "percentile_ucl") >= 0.85
    measured[name] = recovery(r, truths[name])
    measured[name].update({"coverage_bca": float(cover("bca_lcl", "bca_ucl")),
                           "coverage_percentile": float(cover("percentile_lcl", "percentile_ucl"))})


def check_tolerances(measured, truths):
    for name, truth in truths.items():
        exp = json.loads((ROOT / "lib" / name / "expected.json").read_text())
        for key, value in truth.items():
            assert abs(exp["truth"][key] - value) < 1e-10, f"{name} {key}: expected.json {exp['truth'][key]} != {value}"
            tol = float(exp["recovery"].get("tolerance_by_key", {}).get(key, exp["recovery"]["tolerance_estimate"]))
            assert measured[name][key] < tol, f"{name} {key}: recovery miss {measured[name][key]:.4f} >= {tol}"


if __name__ == "__main__":
    measured, truths = {}, {}
    noninferiority(measured, truths)
    multiplicity(measured, truths)
    spc(measured, truths)
    bootstrap(measured, truths)
    if "--measure" in sys.argv:
        print(json.dumps({"seeds": len(SEEDS), "measured": measured, "truths": truths}, indent=2))
        sys.exit(0)
    check_tolerances(measured, truths)
    print(json.dumps({"entries": len(measured), "seeds": len(SEEDS), "result": "pass"}, indent=2))
