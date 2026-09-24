"""Controls and fixed-seed calibration for the instrumental-variable, mediation and case-crossover entries.

instrumental-variables-2sls, mediation-counterfactual-binary, case-crossover-conditional-logistic. Each
fixture regenerates byte for byte and each stated truth is recomputed (mediation by scipy quadrature).
Oracles: the two-stage OLS shortcut reproduces the 2SLS coefficient and gets its standard error wrong;
the mediation sandwich is rebuilt with numerical derivatives; the case-crossover estimate is the exact
conditional MLE found by Newton's method. Negative controls, over 100 seeds: OLS, the difference-method
"indirect effect" and the unmatched logistic regression each miss. Recovery is calibrated by running each
entry's own python.py on 100 generated fixtures. Run with --measure to print the numbers.
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
from scipy import integrate
from scipy.special import expit, logit

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check import parse_harness_block  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
SEEDS = range(int(os.environ.get("MALLARD_SEEDS", "50")))  # 50 in CI; MALLARD_SEEDS=100 reproduces the recorded calibration


def rank_quantile(misses):
    """The recovery statistic: the quantile at the same ORDER STATISTIC as the 99th percentile of 100 draws.

    np.quantile(q=0.99) over 100 misses lands just above the second-largest; over 50 it lands halfway to
    the largest, which is stricter than the calibration the tolerances were set from. The level
    (n - 1.99) / (n - 1) is exactly 0.99 at n = 100 and the same position at any n.
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


def runs(name):
    module = load(name)
    return [run_python(name, list(module.rows(s)))[0] for s in SEEDS]


def miss99(values, truth):
    return rank_quantile([abs(v - truth) for v in values])


# ---------------------------------------------------------------------------------------------------

def iv(measured, truths):
    name = "instrumental-variables-2sls"
    regenerates(name)
    m = load(name)
    truths[name] = {"iv_effect": m.EFFECT}
    out, scope = run_python(name)
    # The two-stage OLS shortcut: same coefficient, wrong residuals, wrong SE.
    Xhat, y = scope["Xhat"], scope["y"]
    b2 = np.linalg.lstsq(Xhat, y, rcond=None)[0]
    assert abs(b2[1] - out["iv_effect"]) < 1e-10
    e2 = y - Xhat @ b2
    se_shortcut = math.sqrt(np.linalg.inv(Xhat.T @ Xhat)[1, 1] * np.sum(e2 ** 2) / (len(y) - 2))
    assert abs(se_shortcut / out["iv_se_conventional"] - 1) > 0.1, "the shortcut SE should be wrong"
    assert out["first_stage_f"] > 100
    r = runs(name)
    assert np.mean([x["ols_effect"] for x in r]) - m.EFFECT > 0.4, "OLS should be confounded upwards"
    measured[name] = {"iv_effect": miss99([x["iv_effect"] for x in r], m.EFFECT),
                      "shortcut_se_ratio": se_shortcut / out["iv_se_conventional"]}


def mediation(measured, truths):
    name = "mediation-counterfactual-binary"
    regenerates(name)
    m = load(name)
    phi = lambda c: math.exp(-0.5 * c * c) / math.sqrt(2 * math.pi)
    risk = lambda a, s: integrate.quad(
        lambda c: (m.p_y(a, 1, c) * m.p_m(s, c) + m.p_y(a, 0, c) * (1 - m.p_m(s, c))) * phi(c),
        -12, 12, epsabs=1e-14, epsrel=1e-13, limit=200)[0]
    r11, r10, r00 = risk(1, 1), risk(1, 0), risk(0, 0)
    truth = {"natural_direct_effect": r10 - r00, "natural_indirect_effect": r11 - r10}
    truths[name] = truth
    stated = m.truth()
    for key, value in truth.items():
        assert abs(stated[key] - value) < 1e-10
    nie_log_or = logit(r11) - logit(r10)

    out, scope = run_python(name)
    # Numerical-derivative sandwich over the same estimating functions.
    Xm, Xy, XM, XY = scope["Xm"], scope["Xy"], scope["XM"], scope["XY"]
    dm, dy = scope["d"]["mediator"].to_numpy(float), scope["d"]["event"].to_numpy(float)
    k1, k2 = Xm.shape[1], Xy.shape[1]
    pairs = scope["pairs"]
    theta = np.r_[scope["g"], scope["b"], scope["mu"]]

    def psi(th):
        g, b, mu = th[:k1], th[k1:k1 + k2], th[k1 + k2:]
        cols = [Xm * (dm - expit(Xm @ g))[:, None], Xy * (dy - expit(Xy @ b))[:, None]]
        for j, (a, s) in enumerate(pairs):
            q = expit(XM(s) @ g)
            cols.append(expit(XY(a, 1) @ b) * q + expit(XY(a, 0) @ b) * (1 - q) - mu[j])
        return np.column_stack(cols)

    A = np.column_stack([(psi(theta + h).sum(0) - psi(theta - h).sum(0)) / 2e-6 for h in np.eye(len(theta)) * 1e-6])
    Ainv = np.linalg.inv(A)
    V = Ainv @ (psi(theta).T @ psi(theta)) @ Ainv.T
    w = np.r_[np.zeros(k1 + k2), -1.0, 0.0, 1.0]
    assert abs(math.sqrt(w @ V @ w) - out["natural_indirect_effect_se"]) < 1e-7
    assert abs(out["natural_direct_effect"] + out["natural_indirect_effect"] - out["total_effect"]) < 1e-12

    r = runs(name)
    diff_method = float(np.mean([x["difference_method_log_odds_ratio"] for x in r]))
    counterfactual = float(np.mean([x["nie_log_odds_ratio"] for x in r]))
    assert abs(counterfactual - nie_log_or) < 0.02, "the counterfactual NIE (log OR) should be unbiased"
    assert abs(diff_method - nie_log_or) > 0.05, "the difference method should miss the indirect effect"
    measured[name] = {key: miss99([x[key] for x in r], value) for key, value in truth.items()}
    measured[name].update({"nie_log_or_truth": nie_log_or, "difference_method_mean": diff_method,
                           "counterfactual_log_or_mean": counterfactual})


def case_crossover(measured, truths):
    name = "case-crossover-conditional-logistic"
    regenerates(name)
    m = load(name)
    truths[name] = {"log_or": m.LOG_OR}
    out, scope = run_python(name)
    d = scope["d"]
    X = d.exposed.to_numpy(float).reshape(-1, m.WINDOWS)
    xc = (X * d.hazard_window.to_numpy().reshape(-1, m.WINDOWS)).sum(1)
    b = 0.0
    for _ in range(100):
        w = np.exp(b * X)
        m1 = (w * X).sum(1) / w.sum(1)
        info = ((w * X * X).sum(1) / w.sum(1) - m1 ** 2).sum()
        b += (xc - m1).sum() / info
    assert abs(b - out["log_or"]) < 1e-8, "not the exact conditional MLE"
    assert abs(1 / math.sqrt(info) - out["log_or_se"]) < 1e-7
    r = runs(name)
    assert np.mean([x["unconditional_log_or"] for x in r]) < m.LOG_OR - 0.1, "ignoring the matching should attenuate"
    measured[name] = {"log_or": miss99([x["log_or"] for x in r], m.LOG_OR)}


def check_tolerances(measured, truths):
    for name, truth in truths.items():
        exp = json.loads((ROOT / "lib" / name / "expected.json").read_text())
        for key, value in truth.items():
            assert abs(exp["truth"][key] - value) < 1e-10, f"{name} {key}: expected.json {exp['truth'][key]} != {value}"
            tol = float(exp["recovery"].get("tolerance_by_key", {}).get(key, exp["recovery"]["tolerance_estimate"]))
            assert measured[name][key] < tol, f"{name} {key}: 99th percentile miss {measured[name][key]:.4f} >= {tol}"


if __name__ == "__main__":
    measured, truths = {}, {}
    iv(measured, truths)
    mediation(measured, truths)
    case_crossover(measured, truths)
    if "--measure" in sys.argv:
        print(json.dumps({"measured": measured, "truths": truths}, indent=2))
        sys.exit(0)
    check_tolerances(measured, truths)
    print(json.dumps({"entries": len(measured), "seeds": len(SEEDS), "result": "pass"}, indent=2))
