"""Controls and fixed-seed calibration for the regression entries in batch C-1.

ancova-baseline-adjusted, binomial-identity-risk-difference, negative-binomial-rate-ratio,
ordinal-proportional-odds. For each: byte-for-byte fixture regeneration; an oracle that maximizes the
model's likelihood directly with scipy (independent of both engines' fitting routines); a negative
control for the pinned trap; every stated truth recomputed from the generator; and recovery measured
over 100 seeds. Run with --measure to print the calibration numbers.
"""
import csv
import importlib.util
import io
import json
import math
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import statsmodels.api as sm
from scipy import optimize, special

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check import parse_harness_block  # noqa: E402

warnings.filterwarnings("ignore", category=sm.tools.sm_exceptions.DomainWarning)
ROOT = Path(__file__).resolve().parent.parent
SEEDS = range(100)


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
    return float(np.quantile([abs(e - truth) for e in estimates], 0.99))


# --- ANCOVA --------------------------------------------------------------------------------------

def ancova_fit(rows):
    y = np.array([r["followup"] for r in rows])
    X = np.column_stack([np.ones(len(rows)), [r["treated"] for r in rows], [r["baseline"] for r in rows]])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    sigma2 = resid @ resid / (len(y) - X.shape[1])
    se = np.sqrt(np.diag(sigma2 * np.linalg.inv(X.T @ X)))
    return beta, se


def ancova(measured):
    name = "ancova-baseline-adjusted"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    rows = list(module.rows())
    beta, se = ancova_fit(rows)
    assert abs(out["effect"] - beta[1]) < 1e-9 and abs(out["effect_se"] - se[1]) < 1e-9
    # Negative controls: the follow-up-only and change-score analyses are less precise.
    t = np.array([r["treated"] for r in rows]); f = np.array([r["followup"] for r in rows])
    b = np.array([r["baseline"] for r in rows])
    def diff_se(v):
        return math.sqrt(v[t == 1].var(ddof=1) / (t == 1).sum() + v[t == 0].var(ddof=1) / (t == 0).sum())
    followup_only, change_score = diff_se(f), diff_se(f - b)
    assert followup_only > 1.2 * out["effect_se"] and change_score > 1.05 * out["effect_se"]
    per_seed = [ancova_fit(list(module.rows(s)))[0] for s in SEEDS]
    measured[name] = {"effect": miss99([p[1] for p in per_seed], module.EFFECT),
                      "baseline_coef": miss99([p[2] for p in per_seed], module.BASELINE_COEF),
                      "se_followup_only": followup_only, "se_change_score": change_score,
                      "se_ancova": out["effect_se"]}


# --- identity-link binomial ----------------------------------------------------------------------

def identity_mle(rows):
    y = np.array([r["outcome"] for r in rows], dtype=float)
    X = np.column_stack([np.ones(len(rows)), [r["treated"] for r in rows], [r["age"] - 60 for r in rows]])
    def nll(b):
        p = np.clip(X @ b, 1e-12, 1 - 1e-12)
        return -(y * np.log(p) + (1 - y) * np.log(1 - p)).sum()
    res = optimize.minimize(nll, np.array([y.mean(), 0.0, 0.0]), method="Nelder-Mead",
                            options={"xatol": 1e-12, "fatol": 1e-14, "maxiter": 20000, "maxfev": 40000})
    return res.x


def identity_rd(measured):
    name = "binomial-identity-risk-difference"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    rows = list(module.rows())
    beta = identity_mle(rows)
    assert abs(out["risk_difference"] - beta[1]) < 1e-5 and abs(out["age_slope"] - beta[2]) < 1e-6
    # Negative control: the default logit link returns a log odds ratio, a different number.
    y = np.array([r["outcome"] for r in rows]); X = sm.add_constant(
        np.column_stack([[r["treated"] for r in rows], [r["age"] - 60 for r in rows]]))
    logit_coef = sm.GLM(y, X, family=sm.families.Binomial()).fit().params[1]
    assert abs(logit_coef - out["risk_difference"]) > 0.1
    per_seed = []
    for s in SEEDS:
        rs = list(module.rows(s))
        yy = np.array([r["outcome"] for r in rs]); XX = sm.add_constant(
            np.column_stack([[r["treated"] for r in rs], [r["age"] - 60 for r in rs]]))
        fit = sm.GLM(yy, XX, family=sm.families.Binomial(link=sm.families.links.Identity())).fit(
            start_params=np.array([yy.mean(), 0.0, 0.0]), tol=1e-12, maxiter=100)
        per_seed.append(fit.params)
    measured[name] = {"risk_difference": miss99([p[1] for p in per_seed], module.RISK_DIFFERENCE),
                      "age_slope": miss99([p[2] for p in per_seed], module.AGE_SLOPE),
                      "logit_coefficient_committed": float(logit_coef)}


# --- negative binomial ---------------------------------------------------------------------------

def nb_mle(rows):
    y = np.array([r["events"] for r in rows], dtype=float)
    t = np.array([r["treated"] for r in rows], dtype=float)
    off = np.log([r["followup_years"] for r in rows])
    def nll(p):
        b0, b1, log_theta = p
        mu = np.exp(b0 + b1 * t + off); th = math.exp(log_theta)
        ll = (special.gammaln(y + th) - special.gammaln(th) - special.gammaln(y + 1)
              + th * np.log(th / (th + mu)) + y * np.log(mu / (th + mu)))
        return -ll.sum()
    res = optimize.minimize(nll, np.array([0.0, 0.0, 0.0]), method="BFGS", options={"gtol": 1e-10})
    return res.x[1], math.exp(res.x[2])


def negative_binomial(measured):
    name = "negative-binomial-rate-ratio"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    rows = list(module.rows())
    b, theta = nb_mle(rows)
    assert abs(out["log_rate_ratio"] - b) < 1e-5 and abs(out["theta"] - theta) < 1e-3
    # Negative controls: a Poisson fit, and a NegativeBinomial family left at its default alpha = 1,
    # give different standard errors.
    y = np.array([r["events"] for r in rows]); X = sm.add_constant(np.array([r["treated"] for r in rows], dtype=float))
    off = np.log([r["followup_years"] for r in rows])
    poisson_se = sm.GLM(y, X, offset=off, family=sm.families.Poisson()).fit().bse[1]
    default_se = sm.GLM(y, X, offset=off, family=sm.families.NegativeBinomial()).fit().bse[1]
    assert poisson_se < 0.85 * out["log_rate_ratio_se"]
    assert abs(default_se - out["log_rate_ratio_se"]) > 1e-3
    per_seed = [nb_mle(list(module.rows(s))) for s in SEEDS]
    measured[name] = {"log_rate_ratio": miss99([p[0] for p in per_seed], module.LOG_RATE_RATIO),
                      "theta": miss99([p[1] for p in per_seed], module.THETA),
                      "poisson_se_committed": float(poisson_se), "default_alpha_se_committed": float(default_se)}


# --- proportional odds ---------------------------------------------------------------------------

def po_mle(rows, levels):
    y = np.array([levels.index(r["severity"]) for r in rows])
    x = np.array([r["treated"] for r in rows], dtype=float)
    k = len(levels)
    def nll(p):
        beta, c1, *incs = p
        cuts = np.concatenate([[c1], c1 + np.cumsum(np.exp(incs))])
        upper = np.append(cuts, np.inf)[y] - beta * x
        lower = np.insert(cuts, 0, -np.inf)[y] - beta * x
        prob = special.expit(upper) - special.expit(lower)
        return -np.log(np.clip(prob, 1e-300, None)).sum()
    res = optimize.minimize(nll, np.zeros(k), method="BFGS", options={"gtol": 1e-10})
    beta, c1, *incs = res.x
    return beta, np.concatenate([[c1], c1 + np.cumsum(np.exp(incs))])


def ordinal(measured):
    name = "ordinal-proportional-odds"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    rows = list(module.rows())
    beta, cuts = po_mle(rows, list(module.LEVELS))
    assert abs(out["log_odds_ratio"] - beta) < 1e-5
    for k, c in enumerate(cuts, 1):
        assert abs(out[f"threshold_{k}"] - c) < 1e-5
    # Negative control: the alphabetical ordering (mild, moderate, none, severe) fits a different model.
    alpha_beta, _ = po_mle(rows, sorted(module.LEVELS))
    assert abs(alpha_beta - out["log_odds_ratio"]) > 0.05
    per_seed = [po_mle(list(module.rows(s)), list(module.LEVELS)) for s in SEEDS]
    truth = {"log_odds_ratio": module.LOG_OR, **{f"threshold_{k}": c for k, c in enumerate(module.CUTS, 1)}}
    measured[name] = {"log_odds_ratio": miss99([p[0] for p in per_seed], module.LOG_OR),
                      **{f"threshold_{k}": miss99([p[1][k - 1] for p in per_seed], truth[f"threshold_{k}"]) for k in (1, 2, 3)},
                      "alphabetical_log_odds_ratio_committed": float(alpha_beta)}


def generator_truth():
    a = load("ancova-baseline-adjusted"); i = load("binomial-identity-risk-difference")
    n = load("negative-binomial-rate-ratio"); o = load("ordinal-proportional-odds")
    return {
        "ancova-baseline-adjusted": {"effect": a.EFFECT, "baseline_coef": a.BASELINE_COEF},
        "binomial-identity-risk-difference": {"risk_difference": i.RISK_DIFFERENCE, "age_slope": i.AGE_SLOPE},
        "negative-binomial-rate-ratio": {"log_rate_ratio": n.LOG_RATE_RATIO, "theta": n.THETA},
        "ordinal-proportional-odds": {"log_odds_ratio": o.LOG_OR, **{f"threshold_{k}": c for k, c in enumerate(o.CUTS, 1)}},
    }


def check_tolerances(measured):
    for name, truth in generator_truth().items():
        stated = json.loads((ROOT / "lib" / name / "expected.json").read_text())["truth"]
        numeric = {k: v for k, v in stated.items() if isinstance(v, (int, float))}
        assert set(numeric) == set(truth), f"{name}: truth keys {sorted(numeric)} != generator {sorted(truth)}"
        for key, value in truth.items():
            assert abs(numeric[key] - value) < 1e-12, f"{name} {key}: expected.json {numeric[key]} != generator {value}"
    for name, misses in measured.items():
        rec = json.loads((ROOT / "lib" / name / "expected.json").read_text())["recovery"]
        for key, miss in misses.items():
            if key not in generator_truth()[name]:
                continue
            tol = float(rec.get("tolerance_by_key", {}).get(key, rec["tolerance_estimate"]))
            assert miss < tol, f"{name} {key}: 99th percentile miss {miss:.4f} >= tolerance {tol}"


if __name__ == "__main__":
    measured = {}
    ancova(measured)
    identity_rd(measured)
    negative_binomial(measured)
    ordinal(measured)
    if "--measure" in sys.argv:
        print(json.dumps(measured, indent=2))
        sys.exit(0)
    check_tolerances(measured)
    print(json.dumps({"entries": len(measured), "seeds": len(SEEDS), "result": "pass"}, indent=2))
