"""Controls and fixed-seed calibration for batch F-1.

propensity-score-matching-att, overlap-weights-ato, g-computation-standardization, aipw-doubly-robust.
Each fixture regenerates byte for byte, and each stated truth is recomputed here with scipy's adaptive
quadrature (the generators use a stdlib Simpson rule). Oracles: the matched pairs are rebuilt by brute
force and checked for the caliper and for no reuse; overlap weights balance the covariate means exactly;
the stacked sandwiches are rebuilt with numerical derivatives; AIPW is shown to be doubly robust by
breaking each working model in turn over 100 seeds. Recovery is calibrated by running each entry's own
python.py on 100 generated fixtures. Run with --measure to print the numbers.
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
import pandas as pd
import statsmodels.api as sm
from scipy import integrate

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check import parse_harness_block  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
SEEDS = range(int(os.environ.get("MALLARD_SEEDS", "50")))  # 50 in CI; MALLARD_SEEDS=100 reproduces the recorded calibration (owner decisions E5, E6)


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
    """Run the entry's python.py in-process, on the committed fixture or on generated rows."""
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


def miss99(estimates, truth):
    return rank_quantile([abs(e - truth) for e in estimates])


def normal_mean(f):
    return integrate.quad(lambda x: f(x) * math.exp(-0.5 * x * x) / math.sqrt(2 * math.pi), -12, 12,
                          epsabs=1e-14, epsrel=1e-13, limit=200)[0]


def calibrate(name, key, truth):
    module = load(name)
    return miss99([run_python(name, list(module.rows(s)))[0][key] for s in SEEDS], truth)


# ---------------------------------------------------------------------------------------------------

def matching(measured, truths):
    name = "propensity-score-matching-att"
    regenerates(name)
    m = load(name)
    num = den = 0.0
    for z, pz in ((0, 1 - m.P_Z), (1, m.P_Z)):
        num += pz * normal_mean(lambda l: m.propensity(l, z) * (m.EFFECT + m.EFFECT_BY_L * l))
        den += pz * normal_mean(lambda l: m.propensity(l, z))
    truths[name] = {"att": num / den}
    assert abs(m.truth_att() - num / den) < 1e-10

    out, scope = run_python(name)
    pairs, dist, caliper = scope["pairs"], scope["dist"], scope["caliper"]
    # Oracle: every pair is inside the caliper, no control is used twice, every pair is treated-control.
    treated_flag = scope["d"]["treated"].to_numpy()
    assert len({c for _, c in pairs}) == len(pairs), "a control was reused"
    assert all(treated_flag[t] == 1 and treated_flag[c] == 0 for t, c in pairs)
    assert max(abs(dist[t] - dist[c]) for t, c in pairs) <= caliper
    # Oracle: each unmatched treated patient had no unused control within the caliper at its turn,
    # checked by replaying the greedy order by brute force.
    used = set()
    matched = dict(pairs)
    controls = np.flatnonzero(treated_flag == 0)
    for t in sorted(np.flatnonzero(treated_flag == 1), key=lambda i: (-dist[i], i)):
        free = [c for c in controls if c not in used]
        best = min(free, key=lambda c: (abs(dist[c] - dist[t]), c))
        if abs(dist[best] - dist[t]) <= caliper:
            assert matched.get(t) == best, "the greedy order was not followed"
            used.add(best)
        else:
            assert t not in matched
    # Balance improves to under 0.1 from over 0.3; matching is not the naive comparison.
    assert out["smd_severity_before"] > 0.3 and abs(out["smd_severity_after"]) < 0.1
    assert out["naive_difference"] - out["att"] > 0.5
    measured[name] = {"att": calibrate(name, "att", num / den), "att_minus_ate": num / den - m.EFFECT}


def overlap(measured, truths):
    name = "overlap-weights-ato"
    regenerates(name)
    m = load(name)
    num = den = 0.0
    for z, pz in ((0, 1 - m.P_Z), (1, m.P_Z)):
        tilt = lambda l, z=z: m.propensity(l, z) * (1 - m.propensity(l, z))
        num += pz * normal_mean(lambda l: tilt(l) * (m.EFFECT + m.EFFECT_BY_L * l))
        den += pz * normal_mean(tilt)
    truths[name] = {"ato": num / den}
    assert abs(m.truth_ato() - num / den) < 1e-10

    out, scope = run_python(name)
    # Exact mean balance: a property of overlap weights with a logistic propensity model.
    assert abs(out["smd_severity_weighted"]) < 1e-10
    # Numerical-derivative sandwich over the same estimating functions.
    X, a, y = scope["X"], scope["a"], scope["y"]
    theta = np.r_[np.asarray(scope["ps"].params), scope["mu1"], scope["mu0"]]
    k = X.shape[1]

    def psi(th):
        e = 1 / (1 + np.exp(-X @ th[:k]))
        w = np.where(a == 1, 1 - e, e)
        return np.column_stack([X * (a - e)[:, None], w * a * (y - th[k]), w * (1 - a) * (y - th[k + 1])])

    A = np.column_stack([(psi(theta + h).sum(0) - psi(theta - h).sum(0)) / 2e-6
                         for h in np.eye(len(theta)) * 1e-6])
    Ainv = np.linalg.inv(A)
    B = psi(theta).T @ psi(theta)
    c = np.r_[np.zeros(k), 1.0, -1.0]
    assert abs(math.sqrt(c @ Ainv @ B @ Ainv.T @ c) - out["ato_se"]) < 1e-7
    assert abs(out["ato_se"] - out["ato_se_weights_fixed"]) > 0.01
    measured[name] = {"ato": calibrate(name, "ato", num / den), "ato_minus_ate": num / den - m.EFFECT}


def gcomp(measured, truths):
    name = "g-computation-standardization"
    regenerates(name)
    m = load(name)
    r1 = r0 = 0.0
    for z, pz in ((0, 1 - m.P_Z), (1, m.P_Z)):
        r1 += pz * normal_mean(lambda l: m.risk(1, l, z))
        r0 += pz * normal_mean(lambda l: m.risk(0, l, z))
    truth = {"risk_difference": r1 - r0, "log_risk_ratio": math.log(r1 / r0)}
    truths[name] = truth
    stated = m.truth()
    for key, value in truth.items():
        assert abs(stated[key] - value) < 1e-10
    # Non-collapsibility in the truth itself: the marginal log OR is below the conditional one.
    marginal = math.log(r1 / (1 - r1)) - math.log(r0 / (1 - r0))
    assert m.OUT[1] - marginal > 0.05

    out, scope = run_python(name)
    X, X1, X0, yv = scope["X"], scope["X1"], scope["X0"], scope["yv"]
    theta = np.r_[np.asarray(scope["fit"].params), scope["mu1"], scope["mu0"]]
    k = X.shape[1]

    def psi(th):
        ex = lambda M: 1 / (1 + np.exp(-M @ th[:k]))
        return np.column_stack([X * (yv - ex(X))[:, None], ex(X1) - th[k], ex(X0) - th[k + 1]])

    A = np.column_stack([(psi(theta + h).sum(0) - psi(theta - h).sum(0)) / 2e-6
                         for h in np.eye(len(theta)) * 1e-6])
    Ainv = np.linalg.inv(A)
    V = Ainv @ (psi(theta).T @ psi(theta)) @ Ainv.T
    c = np.r_[np.zeros(k), 1.0, -1.0]
    assert abs(math.sqrt(c @ V @ c) - out["risk_difference_se"]) < 1e-7
    assert out["log_conditional_odds_ratio"] > out["log_marginal_odds_ratio"]
    assert out["crude_risk_difference"] - out["risk_difference"] > 0.05
    assert out["risk_difference_se_conditional"] < out["risk_difference_se"]
    measured[name] = {key: calibrate(name, key, value) for key, value in truth.items()}


def aipw(measured, truths):
    name = "aipw-doubly-robust"
    regenerates(name)
    m = load(name)
    truth = gcomp_truth = m.truth()["risk_difference"]
    truths[name] = {"ate": truth}
    assert abs(gcomp_truth - truths["g-computation-standardization"]["risk_difference"]) < 1e-12

    # Double robustness, over 100 seeds: break one working model at a time.
    def fit(y, X):
        return sm.Logit(y, X).fit(disp=0, tol=1e-12)

    bad_ps = {"aipw": [], "single": []}
    bad_om = {"aipw": [], "single": []}
    for s in SEEDS:
        d = pd.DataFrame(list(m.rows(s)))
        a, y = d.treated.to_numpy(float), d.event.to_numpy(float)
        one = np.ones(len(d))
        XL = np.column_stack([one, d.severity, d.comorbid])
        XO = np.column_stack([one, a, d.severity, d.comorbid])
        e_good = fit(a, XL).predict(XL)
        e_bad = np.full(len(d), a.mean())  # propensity model with no covariates
        om = fit(y, XO)
        m1 = om.predict(np.column_stack([one, one, d.severity, d.comorbid]))
        m0 = om.predict(np.column_stack([one, 0 * one, d.severity, d.comorbid]))
        omb = fit(y, np.column_stack([one, a]))  # outcome model with no covariates
        mb1 = np.full(len(d), omb.predict(np.array([[1.0, 1.0]]))[0])
        mb0 = np.full(len(d), omb.predict(np.array([[1.0, 0.0]]))[0])

        def est(e, p1, p0):
            return np.mean(p1 + a * (y - p1) / e) - np.mean(p0 + (1 - a) * (y - p0) / (1 - e))

        bad_ps["aipw"].append(est(e_bad, m1, m0))
        bad_ps["single"].append(np.mean(a * y / e_bad) - np.mean((1 - a) * y / (1 - e_bad)))
        bad_om["aipw"].append(est(e_good, mb1, mb0))
        bad_om["single"].append(mb1[0] - mb0[0])
    bias = lambda v: float(np.mean(v) - truth)
    assert abs(bias(bad_ps["aipw"])) < 0.01 and bias(bad_ps["single"]) > 0.08, "propensity model broken"
    assert abs(bias(bad_om["aipw"])) < 0.01 and bias(bad_om["single"]) > 0.08, "outcome model broken"

    out, _ = run_python(name)
    assert abs(out["ate"] - out["gcomp_only"]) < 0.01 and abs(out["ate"] - out["ipw_only"]) < 0.01
    measured[name] = {"ate": calibrate(name, "ate", truth),
                      "bias_bad_ps": {k: bias(v) for k, v in bad_ps.items()},
                      "bias_bad_outcome": {k: bias(v) for k, v in bad_om.items()}}


def check_tolerances(measured, truths):
    for name, truth in truths.items():
        exp = json.loads((ROOT / "lib" / name / "expected.json").read_text())
        for key, value in truth.items():
            assert abs(exp["truth"][key] - value) < 1e-10, f"{name} {key}: expected.json {exp['truth'][key]} != {value}"
            tol = float(exp["recovery"].get("tolerance_by_key", {}).get(key, exp["recovery"]["tolerance_estimate"]))
            assert measured[name][key] < tol, f"{name} {key}: 99th percentile miss {measured[name][key]:.4f} >= {tol}"


if __name__ == "__main__":
    measured, truths = {}, {}
    matching(measured, truths)
    overlap(measured, truths)
    gcomp(measured, truths)
    aipw(measured, truths)
    if "--measure" in sys.argv:
        print(json.dumps({"measured": measured, "truths": truths}, indent=2))
        sys.exit(0)
    check_tolerances(measured, truths)
    print(json.dumps({"entries": len(measured), "seeds": len(SEEDS), "result": "pass"}, indent=2))
