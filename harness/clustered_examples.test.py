"""Controls and fixed-seed calibration for batch E-2.

stepped-wedge-mixed-model, cluster-robust-cr2-small-sample (R and Python), mixed-effects-poisson-clustered
(R only). The stepped-wedge oracle is an independent REML fit in numpy (profiling the variance ratio);
CR2 is checked against R's clubSandwich by the harness and against CR1 and the naive SE here; the mixed
Poisson is checked in R against adaptive quadrature, and its recovery is measured by fitting 100 fixtures
in one Rscript call. Run with --measure to print the numbers. Locally, MALLARD_R_LIBS may name an R
library to put first on the path.
"""
import csv
import importlib.util
import io
import json
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check import parse_harness_block  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
SEEDS = range(100)


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "lib" / name / "fixture.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_csv(path, generated):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(generated[0]))
        writer.writeheader()
        writer.writerows(generated)


def regenerates(name):
    generated = list(load(name).rows())
    out = io.StringIO(newline="")
    writer = csv.DictWriter(out, fieldnames=list(generated[0]))
    writer.writeheader()
    writer.writerows(generated)
    assert out.getvalue().encode() == (ROOT / "lib" / name / "fixture.csv").read_bytes(), \
        f"{name}: fixture.csv does not regenerate byte for byte"


def run_r(code, cwd=None):
    libs = os.environ.get("MALLARD_R_LIBS")
    prelude = f'.libPaths(c("{libs}", .libPaths())); ' if libs else ""
    return subprocess.check_output(["Rscript", "-e", prelude + code], cwd=cwd, text=True, stderr=subprocess.DEVNULL)


def python_output(name):
    return parse_harness_block(subprocess.check_output(
        [sys.executable, "python.py"], cwd=ROOT / "lib" / name, text=True))


def miss99(estimates, truth):
    return float(np.quantile([abs(e - truth) for e in estimates], 0.99))


# --- stepped wedge: an independent REML ----------------------------------------------------------

def reml_random_intercept(X, y, groups):
    """REML for y = X b + u_group + e, profiling sigma_e^2 and maximizing over lambda = s2u / s2e."""
    blocks = [np.flatnonzero(groups == g) for g in np.unique(groups)]
    n, p = X.shape

    def pieces(lam):
        xtvx = np.zeros((p, p)); xtvy = np.zeros(p); logdet = 0.0; parts = []
        for idx in blocks:
            m = len(idx)
            Vinv = np.eye(m) - lam / (1 + lam * m) * np.ones((m, m))   # (I + lam J)^-1
            xtvx += X[idx].T @ Vinv @ X[idx]; xtvy += X[idx].T @ Vinv @ y[idx]
            logdet += math.log(1 + lam * m)
            parts.append((idx, Vinv))
        b = np.linalg.solve(xtvx, xtvy)
        rss = sum((y[idx] - X[idx] @ b) @ Vinv @ (y[idx] - X[idx] @ b) for idx, Vinv in parts)
        return b, xtvx, logdet, rss

    def neg(loglam):
        b, xtvx, logdet, rss = pieces(math.exp(loglam))
        s2 = rss / (n - p)
        return 0.5 * ((n - p) * math.log(s2) + logdet + np.linalg.slogdet(xtvx)[1])

    res = optimize.minimize_scalar(neg, bounds=(-12, 6), method="bounded", options={"xatol": 1e-12})
    lam = math.exp(res.x)
    b, xtvx, logdet, rss = pieces(lam)
    s2e = rss / (n - p)
    return b, np.sqrt(np.diag(np.linalg.inv(xtvx)) * s2e), math.sqrt(lam * s2e), math.sqrt(s2e)


def sw_design(d, with_period=True):
    cols = [np.ones(len(d)), d.intervention.to_numpy(float)]
    if with_period:
        cols += [(d.period == p).to_numpy(float) for p in sorted(d.period.unique())[1:]]
    return np.column_stack(cols)


def stepped_wedge(measured):
    name = "stepped-wedge-mixed-model"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    d = pd.DataFrame(list(module.rows()))
    b, se, sd_ward, sd_resid = reml_random_intercept(sw_design(d), d.y.to_numpy(), d.ward.to_numpy())
    assert abs(out["effect"] - b[1]) < 1e-6 and abs(out["effect_se"] - se[1]) < 1e-6
    assert abs(out["sd_ward"] - sd_ward) < 1e-5 and abs(out["sd_residual"] - sd_resid) < 1e-6
    # Negative control: without period effects the secular trend is credited to the intervention.
    b0, *_ = reml_random_intercept(sw_design(d, with_period=False), d.y.to_numpy(), d.ward.to_numpy())
    assert abs(b0[1] - out["effect"]) > 0.5
    per_seed = []
    for s in SEEDS:
        ds = pd.DataFrame(list(module.rows(s)))
        bs, _, sw, sr = reml_random_intercept(sw_design(ds), ds.y.to_numpy(), ds.ward.to_numpy())
        per_seed.append((bs[1], sw, sr))
    measured[name] = {"effect": miss99([p[0] for p in per_seed], module.EFFECT),
                      "sd_ward": miss99([p[1] for p in per_seed], 2.0),
                      "sd_residual": miss99([p[2] for p in per_seed], 8.0),
                      "no_period_effect_committed": float(b0[1])}


# --- CR2 -----------------------------------------------------------------------------------------

def cr2(measured):
    name = "cluster-robust-cr2-small-sample"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    d = pd.DataFrame(list(module.rows()))
    X = np.column_stack([np.ones(len(d)), d.treated, d.age - 50]).astype(float)
    y = d.y.to_numpy()
    import statsmodels.api as sm
    ols = sm.OLS(y, X).fit()
    assert abs(out["effect"] - ols.params[1]) < 1e-9
    # Negative controls: the CR1-style cluster sandwich (statsmodels' default cluster option) and the
    # naive OLS SE are both smaller than CR2, and CR2's Satterthwaite df sit below G - 1.
    cr1 = sm.OLS(y, X).fit(cov_type="cluster", cov_kwds={"groups": d.clinic.to_numpy()}).bse[1]
    assert cr1 < out["effect_cr2_se"] and ols.bse[1] < out["effect_cr2_se"]
    assert out["effect_df"] < d.clinic.nunique() - 1
    est = []
    for s in SEEDS:
        ds = pd.DataFrame(list(module.rows(s)))
        Xs = np.column_stack([np.ones(len(ds)), ds.treated, ds.age - 50]).astype(float)
        est.append(np.linalg.lstsq(Xs, ds.y.to_numpy(), rcond=None)[0][1])
    measured[name] = {"effect": miss99(est, module.EFFECT), "cr1_se_committed": float(cr1),
                      "naive_se_committed": float(ols.bse[1])}


# --- mixed Poisson (R only) ----------------------------------------------------------------------

def mixed_poisson(measured):
    name = "mixed-effects-poisson-clustered"
    regenerates(name)
    module = load(name)
    here = ROOT / "lib" / name
    out = parse_harness_block(run_r('source("r.R")', cwd=here))
    # Oracle: adaptive quadrature (nAGQ = 10) against the Laplace fit the entry reports.
    agq = run_r('suppressPackageStartupMessages(library(lme4)); d <- read.csv("fixture.csv"); '
                'f <- glmer(infections ~ intervention + offset(log(days)) + (1 | ward), data = d, family = poisson, nAGQ = 10); '
                'g <- glm(infections ~ intervention + offset(log(days)), data = d, family = poisson); '
                'cat(sprintf("%.12f %.12f\\n", fixef(f)[2], sqrt(vcov(g)[2, 2])))', cwd=here).split()
    assert abs(float(agq[0]) - out["log_rate_ratio"]) < 1e-3
    # Negative control: ignoring the wards gives a smaller standard error.
    assert float(agq[1]) < 0.95 * out["log_rate_ratio_se"]
    with tempfile.TemporaryDirectory() as tmp:
        for s in SEEDS:
            write_csv(Path(tmp) / f"seed_{s}.csv", list(module.rows(s)))
        est = [float(v) for v in run_r(
            f'suppressPackageStartupMessages(library(lme4)); for (s in 0:{len(SEEDS) - 1}) {{ '
            f'd <- read.csv(file.path("{tmp}", paste0("seed_", s, ".csv"))); '
            f'f <- glmer(infections ~ intervention + offset(log(days)) + (1 | ward), data = d, family = poisson); '
            f'cat(sprintf("%.12f\\n", fixef(f)[2])) }}').split()]
    measured[name] = {"log_rate_ratio": miss99(est, module.LOG_RATE_RATIO),
                      "agq10_minus_laplace": float(agq[0]) - out["log_rate_ratio"],
                      "glm_se_committed": float(agq[1])}


def generator_truth():
    return {
        "stepped-wedge-mixed-model": {"effect": load("stepped-wedge-mixed-model").EFFECT},
        "cluster-robust-cr2-small-sample": {"effect": load("cluster-robust-cr2-small-sample").EFFECT},
        "mixed-effects-poisson-clustered": {"log_rate_ratio": load("mixed-effects-poisson-clustered").LOG_RATE_RATIO},
    }


def check_tolerances(measured):
    truths = generator_truth()
    for name, truth in truths.items():
        stated = json.loads((ROOT / "lib" / name / "expected.json").read_text())["truth"]
        numeric = {k: v for k, v in stated.items() if isinstance(v, (int, float))}
        for key, value in truth.items():
            assert abs(numeric[key] - value) < 1e-12, f"{name} {key}: expected.json {numeric[key]} != generator {value}"
    for name, misses in measured.items():
        rec = json.loads((ROOT / "lib" / name / "expected.json").read_text())["recovery"]
        stated = json.loads((ROOT / "lib" / name / "expected.json").read_text())["truth"]
        for key, miss in misses.items():
            if key not in stated:
                continue
            tol = float(rec.get("tolerance_by_key", {}).get(key, rec["tolerance_estimate"]))
            assert miss < tol, f"{name} {key}: 99th percentile miss {miss:.4f} >= tolerance {tol}"


if __name__ == "__main__":
    measured = {}
    stepped_wedge(measured)
    cr2(measured)
    mixed_poisson(measured)
    if "--measure" in sys.argv:
        print(json.dumps(measured, indent=2))
        sys.exit(0)
    check_tolerances(measured)
    print(json.dumps({"entries": len(measured), "seeds": len(SEEDS), "result": "pass"}, indent=2))
