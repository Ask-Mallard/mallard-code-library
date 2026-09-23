"""Controls and fixed-seed calibration for batch C-2.

multinomial-logistic-regression, quantile-regression-median, restricted-cubic-spline,
zero-inflated-poisson. For each: byte-for-byte fixture regeneration; an oracle computed here by a
route independent of the entry's Python; a negative control for the pinned trap; every stated truth
recomputed from the generator; recovery measured over 100 seeds. Run with --measure to print numbers.
"""
import csv
import importlib.util
import io
import json
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import statsmodels.api as sm
from scipy import optimize, special
from scipy.optimize import linprog
from statsmodels.discrete.count_model import ZeroInflatedPoisson
from statsmodels.tools.numdiff import approx_hess3

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


def regenerates(name):
    generated = list(load(name).rows())
    out = io.StringIO(newline="")
    writer = csv.DictWriter(out, fieldnames=list(generated[0]))
    writer.writeheader()
    writer.writerows(generated)
    assert out.getvalue().encode() == (ROOT / "lib" / name / "fixture.csv").read_bytes(), \
        f"{name}: fixture.csv does not regenerate byte for byte"


def python_output(name):
    return parse_harness_block(subprocess.check_output(
        [sys.executable, "python.py"], cwd=ROOT / "lib" / name, text=True, stderr=subprocess.DEVNULL))


def miss99(estimates, truth):
    return rank_quantile([abs(e - truth) for e in estimates])


# --- multinomial ---------------------------------------------------------------------------------

def mn_mle(rows, levels):
    y = np.array([levels.index(r["destination"]) for r in rows])
    t = np.array([r["treated"] for r in rows], dtype=float)
    def nll(p):
        eta = np.column_stack([np.zeros(len(y))] + [p[2 * j] + p[2 * j + 1] * t for j in range(len(levels) - 1)])
        return -(eta[np.arange(len(y)), y] - special.logsumexp(eta, axis=1)).sum()
    res = optimize.minimize(nll, np.zeros(2 * (len(levels) - 1)), method="BFGS", options={"gtol": 1e-10})
    return res.x


def multinomial(measured):
    name = "multinomial-logistic-regression"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    rows = list(module.rows())
    levels = ["home", "rehab", "nursing"]
    p = mn_mle(rows, levels)
    assert abs(out["rehab_log_rrr"] - p[1]) < 1e-5 and abs(out["nursing_log_rrr"] - p[3]) < 1e-5
    # Negative control: with rehab as the reference, the "nursing" coefficient answers another question.
    alt = mn_mle(rows, ["rehab", "home", "nursing"])
    assert abs(alt[3] - out["nursing_log_rrr"]) > 0.1
    per_seed = [mn_mle(list(module.rows(s)), levels) for s in SEEDS]
    measured[name] = {"rehab_log_rrr": miss99([q[1] for q in per_seed], module.COEF["rehab"][1]),
                      "nursing_log_rrr": miss99([q[3] for q in per_seed], module.COEF["nursing"][1]),
                      "rehab_intercept": miss99([q[0] for q in per_seed], module.COEF["rehab"][0]),
                      "nursing_intercept": miss99([q[2] for q in per_seed], module.COEF["nursing"][0]),
                      "nursing_vs_rehab_reference_committed": float(alt[3])}


# --- quantile regression -------------------------------------------------------------------------

def median_fit(rows):
    X = np.column_stack([np.ones(len(rows)), [r["treated"] for r in rows], [r["age"] - 60 for r in rows]])
    y = np.array([r["los_days"] for r in rows])
    n, p = X.shape
    res = linprog(np.concatenate([np.zeros(p), np.full(2 * n, 0.5)]), A_eq=np.hstack([X, np.eye(n), -np.eye(n)]),
                  b_eq=y, bounds=[(None, None)] * p + [(0, None)] * (2 * n), method="highs")
    return res.x[:p], X, y


def quantile(measured):
    name = "quantile-regression-median"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    rows = list(module.rows())
    beta, X, y = median_fit(rows)
    # Oracle by a different algorithm: statsmodels' iteratively reweighted QuantReg, to its precision.
    qr = sm.QuantReg(y, X).fit(q=0.5, max_iter=5000, p_tol=1e-10)
    assert abs(qr.params[1] - out["median_effect"]) < 1e-4
    # Optimality: the committed estimate's check loss is no worse than nearby points.
    loss = lambda b: np.abs(y - X @ b).sum()
    for step in (1e-3, -1e-3):
        assert loss(beta) <= loss(beta + np.array([0, step, 0])) + 1e-9
    # Negative control: least squares estimates the MEAN effect (true 1.80), a different number.
    ols = sm.OLS(y, X).fit().params[1]
    assert abs(ols - out["median_effect"]) > 0.1
    per_seed = [median_fit(list(module.rows(s)))[0] for s in SEEDS]
    measured[name] = {"median_effect": miss99([b[1] for b in per_seed], module.MEDIAN_EFFECT),
                      "age_slope": miss99([b[2] for b in per_seed], module.AGE_SLOPE),
                      "ols_effect_committed": float(ols)}


# --- restricted cubic spline ---------------------------------------------------------------------

def spline_predictions(rows, knots):
    def basis(x):
        x = np.asarray(x, dtype=float)
        t = knots; k = len(t)
        pos = lambda v: np.maximum(v, 0.0) ** 3
        cols = [np.ones_like(x), x]
        # Harrell's normalization by (t_k - t_1)^2, unlike python.py: a different scaling of one basis.
        for j in range(k - 2):
            cols.append((pos(x - t[j]) - pos(x - t[k - 2]) * (t[k - 1] - t[j]) / (t[k - 1] - t[k - 2])
                         + pos(x - t[k - 1]) * (t[k - 2] - t[j]) / (t[k - 1] - t[k - 2])) / (t[-1] - t[0]) ** 2)
        return np.column_stack(cols)
    x = np.array([r["age"] for r in rows]); y = np.array([r["sbp"] for r in rows])
    b, *_ = np.linalg.lstsq(basis(x), y, rcond=None)
    return basis([40.0, 50.0, 70.0]) @ b


def spline(measured):
    name = "restricted-cubic-spline"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    rows = list(module.rows())
    pred = spline_predictions(rows, np.array(module.KNOTS))
    assert abs(out["pred_40"] - pred[0]) < 1e-8 and abs(out["contrast_70_50"] - (pred[2] - pred[1])) < 1e-8
    # Negative control: knots at the data's quantiles (what ns(age, df = 4) does) fit another curve.
    ages = np.array([r["age"] for r in rows])
    q_knots = np.quantile(ages, [0.05, 0.275, 0.5, 0.725, 0.95])
    q_pred = spline_predictions(rows, q_knots)
    assert abs((q_pred[2] - q_pred[1]) - out["contrast_70_50"]) > 1e-3
    truth = {"pred_40": module.curve(40), "pred_50": module.curve(50), "pred_70": module.curve(70),
             "contrast_70_50": module.curve(70) - module.curve(50)}
    per_seed = [spline_predictions(list(module.rows(s)), np.array(module.KNOTS)) for s in SEEDS]
    measured[name] = {"pred_40": miss99([p[0] for p in per_seed], truth["pred_40"]),
                      "pred_50": miss99([p[1] for p in per_seed], truth["pred_50"]),
                      "pred_70": miss99([p[2] for p in per_seed], truth["pred_70"]),
                      "contrast_70_50": miss99([p[2] - p[1] for p in per_seed], truth["contrast_70_50"]),
                      "quantile_knot_contrast_committed": float(q_pred[2] - q_pred[1])}


# --- zero-inflated Poisson -----------------------------------------------------------------------

def zip_loglik(q, X, y):
    pi = special.expit(q[0]); mu = np.exp(X @ q[1:])
    l0 = np.log(pi + (1 - pi) * np.exp(-mu))
    lp = np.log(1 - pi) - mu + y * np.log(mu) - special.gammaln(y + 1)
    return np.where(y == 0, l0, lp).sum()


def zip_mle(rows):
    X = np.column_stack([np.ones(len(rows)), [r["treated"] for r in rows]]).astype(float)
    y = np.array([r["visits"] for r in rows], dtype=float)
    res = optimize.minimize(lambda q: -zip_loglik(q, X, y), np.zeros(3), method="BFGS", options={"gtol": 1e-10})
    return res.x, X, y


def zero_inflated(measured):
    name = "zero-inflated-poisson"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    rows = list(module.rows())
    q, X, y = zip_mle(rows)
    assert abs(out["log_rate_ratio"] - q[2]) < 1e-5 and abs(out["inflation_logit"] - q[0]) < 1e-5
    se = np.sqrt(np.diag(np.linalg.inv(-approx_hess3(q, lambda p: zip_loglik(p, X, y)))))
    assert abs(out["log_rate_ratio_se"] - se[2]) < 1e-5
    # Negative control, and a standing record: statsmodels' own bse for this model is smaller than the
    # log-likelihood's curvature. If a future statsmodels changes this, the control says so.
    model = ZeroInflatedPoisson(y, X, exog_infl=np.ones((len(y), 1)), inflation="logit")
    fit = model.fit(method="newton", maxiter=200, tol=1e-12, disp=False)
    assert abs(fit.bse[2] - out["log_rate_ratio_se"]) > 5e-4
    per_seed = [zip_mle(list(module.rows(s)))[0] for s in SEEDS]
    measured[name] = {"log_rate_ratio": miss99([p[2] for p in per_seed], module.LOG_RATE_RATIO),
                      "inflation_logit": miss99([p[0] for p in per_seed], module.INFLATION_LOGIT),
                      "count_intercept": miss99([p[1] for p in per_seed], module.LOG_MEAN_CONTROL),
                      "statsmodels_bse_committed": float(fit.bse[2])}


def generator_truth():
    m = load("multinomial-logistic-regression"); q = load("quantile-regression-median")
    s = load("restricted-cubic-spline"); z = load("zero-inflated-poisson")
    return {
        "multinomial-logistic-regression": {"rehab_log_rrr": m.COEF["rehab"][1], "nursing_log_rrr": m.COEF["nursing"][1],
                                            "rehab_intercept": m.COEF["rehab"][0], "nursing_intercept": m.COEF["nursing"][0]},
        "quantile-regression-median": {"median_effect": q.MEDIAN_EFFECT, "age_slope": q.AGE_SLOPE},
        "restricted-cubic-spline": {"pred_40": s.curve(40), "pred_50": s.curve(50), "pred_70": s.curve(70),
                                    "contrast_70_50": s.curve(70) - s.curve(50)},
        "zero-inflated-poisson": {"log_rate_ratio": z.LOG_RATE_RATIO, "inflation_logit": z.INFLATION_LOGIT,
                                  "count_intercept": z.LOG_MEAN_CONTROL},
    }


def check_tolerances(measured):
    truths = generator_truth()
    for name, truth in truths.items():
        stated = json.loads((ROOT / "lib" / name / "expected.json").read_text())["truth"]
        numeric = {k: v for k, v in stated.items() if isinstance(v, (int, float))}
        assert set(numeric) == set(truth), f"{name}: truth keys {sorted(numeric)} != generator {sorted(truth)}"
        for key, value in truth.items():
            assert abs(numeric[key] - value) < 1e-12, f"{name} {key}: expected.json {numeric[key]} != generator {value}"
    for name, misses in measured.items():
        rec = json.loads((ROOT / "lib" / name / "expected.json").read_text())["recovery"]
        for key, miss in misses.items():
            if key not in truths[name]:
                continue
            tol = float(rec.get("tolerance_by_key", {}).get(key, rec["tolerance_estimate"]))
            assert miss < tol, f"{name} {key}: 99th percentile miss {miss:.4f} >= tolerance {tol}"


if __name__ == "__main__":
    measured = {}
    multinomial(measured)
    quantile(measured)
    spline(measured)
    zero_inflated(measured)
    if "--measure" in sys.argv:
        print(json.dumps(measured, indent=2))
        sys.exit(0)
    check_tolerances(measured)
    print(json.dumps({"entries": len(measured), "seeds": len(SEEDS), "result": "pass"}, indent=2))
