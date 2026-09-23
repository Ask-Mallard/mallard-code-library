"""Controls and fixed-seed calibration for batch H-2 (agreement, reliability and decision curves).

weighted-kappa-agreement, intraclass-correlation, bland-altman-limits, decision-curve-analysis. Each
fixture regenerates byte for byte and each truth is recomputed independently (kappa from the joint
distribution as a matrix product; net benefit with scipy quadrature). Oracles: the ICC mean squares are
rebuilt from explicit sums of squares; the limits of agreement from the mean and SD; net benefit from a
confusion table. Negative controls over the seeds: percent agreement exceeds kappa; consistency exceeds
absolute agreement when raters differ; the methods correlate above 0.9 while their bias is clearly not
0; the model beats treating everyone at thresholds 0.2 and 0.3.

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
from scipy import integrate, stats
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

def kappa(measured, truths):
    name = "weighted-kappa-agreement"
    regenerates(name)
    m = load(name)
    K = m.K
    report = np.array([[m.p_report(r, t) for t in range(1, K + 1)] for r in range(1, K + 1)])  # [r, t]
    joint = report @ np.diag(m.P_TRUE) @ report.T
    lev = np.arange(K)
    truth = {}
    for key, w in (("kappa", np.eye(K)), ("kappa_quadratic", 1 - (lev[:, None] - lev[None, :]) ** 2 / (K - 1) ** 2)):
        po, pe = np.sum(w * joint), np.sum(w * np.outer(joint.sum(1), joint.sum(0)))
        truth[key] = (po - pe) / (1 - pe)
    stated = m.truth()
    for key, value in truth.items():
        assert abs(stated[key] - value) < 1e-12
    truths[name] = truth
    r = runs(name)
    assert all(x["percent_agreement"] > x["kappa"] for x in r), "percent agreement should exceed kappa"
    assert np.mean([x["kappa_quadratic"] - x["kappa"] for x in r]) > 0.2, "weights should matter on this scale"
    measured[name] = recovery(r, truth)


def icc(measured, truths):
    name = "intraclass-correlation"
    regenerates(name)
    m = load(name)
    truths[name] = m.truth()
    assert abs(truths[name]["icc_agreement"] - 64 / 89) < 1e-12 and abs(truths[name]["icc_consistency"] - 0.8) < 1e-12
    out, s = run_python(name)
    Y = s["Y"]
    ns, nr = Y.shape
    grand = Y.mean()
    ss_rows = nr * np.sum((Y.mean(axis=1) - grand) ** 2)
    ss_cols = ns * np.sum((Y.mean(axis=0) - grand) ** 2)
    ss_err = np.sum((Y - Y.mean(axis=1, keepdims=True) - Y.mean(axis=0, keepdims=True) + grand) ** 2)
    msr, msc, mse = ss_rows / (ns - 1), ss_cols / (nr - 1), ss_err / ((ns - 1) * (nr - 1))
    icc_a = (msr - mse) / (msr + (nr - 1) * mse + nr / ns * (msc - mse))
    assert abs(icc_a - out["icc_agreement"]) < 1e-9 and abs((msr - mse) / (msr + (nr - 1) * mse) - out["icc_consistency"]) < 1e-9
    r = runs(name)
    assert np.mean([x["icc_consistency"] - x["icc_agreement"] for x in r]) > 0.03, "rater offsets should cost agreement"
    measured[name] = recovery(r, truths[name])


def bland_altman(measured, truths):
    name = "bland-altman-limits"
    regenerates(name)
    m = load(name)
    t = m.truth()
    truths[name] = {"bias": t["bias"], "loa_lower": t["loa_lower"], "loa_upper": t["loa_upper"]}
    out, s = run_python(name)
    dif = s["dif"]
    z = stats.norm.ppf(0.975)
    assert abs(dif.mean() - z * dif.std(ddof=1) - out["loa_lower"]) < 1e-9
    r = runs(name)
    assert all(x["correlation"] > 0.9 for x in r) and np.mean([x["bias_lcl"] > 0 for x in r]) > 0.9, \
        "high correlation should coexist with a clear bias"
    measured[name] = recovery(r, truths[name])


def decision_curve(measured, truths):
    name = "decision-curve-analysis"
    regenerates(name)
    m = load(name)
    phi = stats.norm.pdf
    risk = lambda x: expit(m.INTERCEPT + m.SLOPE * x)
    prev = integrate.quad(lambda x: risk(x) * phi(x), -12, 12, epsabs=1e-13)[0]
    truth = {}
    for thr in m.THRESHOLDS:
        cut = (math.log(thr / (1 - thr)) - m.INTERCEPT) / m.SLOPE
        tp = integrate.quad(lambda x: risk(x) * phi(x), cut, 12, epsabs=1e-13)[0]
        fp = integrate.quad(lambda x: (1 - risk(x)) * phi(x), cut, 12, epsabs=1e-13)[0]
        key = f"{thr:.1f}".replace(".", "")
        truth[f"net_benefit_model_t{key}"] = tp - fp * thr / (1 - thr)
        truth[f"net_benefit_all_t{key}"] = prev - (1 - prev) * thr / (1 - thr)
    stated = m.truth()
    for key, value in truth.items():
        assert abs(stated[key] - value) < 1e-8, (key, stated[key], value)
    truths[name] = truth
    out, s = run_python(name)
    d = s["d"]
    treat = (d.predicted_risk > 0.2).to_numpy()
    ev = d.event.to_numpy()
    table = np.array([[np.sum(treat & (ev == 1)), np.sum(treat & (ev == 0))]])
    assert abs((table[0, 0] - table[0, 1] * 0.25) / len(d) - out["net_benefit_model_t02"]) < 1e-9
    r = runs(name)
    for key in ("02", "03"):
        assert all(x[f"net_benefit_model_t{key}"] > x[f"net_benefit_all_t{key}"] for x in r), "model should beat treat-all"
    measured[name] = recovery(r, truth)


def check_tolerances(measured, truths):
    for name, truth in truths.items():
        exp = json.loads((ROOT / "lib" / name / "expected.json").read_text())
        for key, value in truth.items():
            assert abs(exp["truth"][key] - value) < 1e-8, f"{name} {key}: expected.json {exp['truth'][key]} != {value}"
            tol = float(exp["recovery"].get("tolerance_by_key", {}).get(key, exp["recovery"]["tolerance_estimate"]))
            assert measured[name][key] < tol, f"{name} {key}: recovery miss {measured[name][key]:.4f} >= {tol}"


if __name__ == "__main__":
    measured, truths = {}, {}
    kappa(measured, truths)
    icc(measured, truths)
    bland_altman(measured, truths)
    decision_curve(measured, truths)
    if "--measure" in sys.argv:
        print(json.dumps({"seeds": len(SEEDS), "measured": measured, "truths": truths}, indent=2))
        sys.exit(0)
    check_tolerances(measured, truths)
    print(json.dumps({"entries": len(measured), "seeds": len(SEEDS), "result": "pass"}, indent=2))
