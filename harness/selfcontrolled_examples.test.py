"""Controls and fixed-seed calibration for batch F-4.

self-controlled-case-series, nested-case-control-risk-set (R and Python), target-trial-clone-censor-weight
(R only, truth check only by owner decision E3). Each fixture regenerates byte for byte.

Oracles: the SCCS estimate is refitted from the conditional (multinomial) likelihood directly; the nested
case-control estimate is the exact conditional MLE by Newton's method; the target-trial analysis is
re-implemented in numpy and must reproduce R's committed output, and its truth is recomputed by a second
Monte Carlo run with a different seed. Negative controls over 100 seeds: SCCS without age, controls drawn
by cumulative sampling when the outcome is common, and the unweighted and ever-versus-never comparisons.
Recovery: each entry's python.py (or, for the R-only entry, r.R in one Rscript call) on 100 fixtures.
Run with --measure to print the numbers. Locally, MALLARD_R_LIBS may name an R library to put first.
"""
import contextlib
import csv
import importlib.util
import io
import json
import math
import os
import random
import runpy
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import optimize
from scipy.special import expit

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check import parse_harness_block  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
SEEDS = range(100)


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
    """Run the entry's r.R once per generated fixture, in one Rscript call; returns parsed outputs."""
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
    return float(np.quantile([abs(v - truth) for v in values], 0.99))


# ---------------------------------------------------------------------------------------------------

def sccs(measured, truths):
    name = "self-controlled-case-series"
    regenerates(name)
    m = load(name)
    truths[name] = {"log_irr": m.LOG_RR}
    out, scope = run_python(name)
    d = scope["d"]
    # The conditional likelihood: each child's events multinomial over their intervals.
    risk = d.risk.to_numpy(float)
    band2, band3 = (d.age_band == 2).to_numpy(float), (d.age_band == 3).to_numpy(float)
    logdays, events, child = np.log(d.days.to_numpy(float)), d.events.to_numpy(float), d.child.to_numpy()
    starts = np.r_[0, np.flatnonzero(np.diff(child)) + 1]

    def negll(th):
        eta = logdays + th[0] * risk + th[1] * band2 + th[2] * band3
        lse = np.logaddexp.reduceat(eta, starts)
        return -(np.sum(events * eta) - np.sum(np.add.reduceat(events, starts) * lse))

    fit = optimize.minimize(negll, np.zeros(3), method="BFGS", options={"gtol": 1e-10})
    assert abs(fit.x[0] - out["log_irr"]) < 1e-6, "Poisson with child effects is not the conditional SCCS fit"
    runs = [run_python(name, list(m.rows(s)))[0] for s in SEEDS]
    assert np.mean([r["log_irr_age_ignored"] for r in runs]) - m.LOG_RR > 0.2, "ignoring age should inflate it"
    measured[name] = {"log_irr": miss99([r["log_irr"] for r in runs], m.LOG_RR)}


def ncc(measured, truths):
    name = "nested-case-control-risk-set"
    regenerates(name)
    m = load(name)
    truths[name] = {"log_hr": m.LOG_HR}
    out, scope = run_python(name)
    d = scope["d"].sort_values(["set", "case"], ascending=[True, False])
    X = d[["exposed", "age"]].to_numpy(float).reshape(-1, 1 + m.CONTROLS, 2)
    beta = np.zeros(2)
    for _ in range(100):
        w = np.exp(X @ beta)
        p = w / w.sum(1, keepdims=True)
        mean = (p[..., None] * X).sum(1)
        grad = (X[:, 0, :] - mean).sum(0)
        info = sum(np.einsum("j,jk,jl->kl", p[s], X[s], X[s]) - np.outer(mean[s], mean[s]) for s in range(len(X)))
        beta += np.linalg.solve(info, grad)
    assert abs(beta[0] - out["log_hr"]) < 1e-8, "not the exact conditional MLE"
    runs = [run_python(name, list(m.rows(s)))[0] for s in SEEDS]
    # Cumulative sampling from the same cohorts: every member event-free at the end of follow-up as a
    # control (with a common outcome there are not four per case to sample).
    cumulative = []
    for s in SEEDS:
        people = m.cohort(random.Random(s))
        cases = [p for p in people if p[1]]
        controls = [p for p in people if not p[1]]
        rows = cases + controls
        y = np.r_[np.ones(len(cases)), np.zeros(len(controls))]
        Xc = sm.add_constant(np.array([[p[2], p[3]] for p in rows], float))
        cumulative.append(sm.Logit(y, Xc).fit(disp=0).params[1])
    assert abs(np.mean([r["log_hr"] for r in runs]) - m.LOG_HR) < 0.03, "risk-set sampling should be unbiased"
    assert np.mean(cumulative) - m.LOG_HR > 0.15, "cumulative sampling should not estimate the hazard ratio"
    measured[name] = {"log_hr": miss99([r["log_hr"] for r in runs], m.LOG_HR),
                      "cumulative_sampling_mean_log_or": float(np.mean(cumulative))}


def ccw_numpy(d, grace=3, intervals=10):
    d = d.sort_values(["id", "interval"]).reset_index(drop=True)
    s = d.assign(t=np.where(d.treated == 1, d.interval, np.inf)).groupby("id")["t"].transform("min").to_numpy()
    k = d.interval.to_numpy()
    elig = (k < grace) & (s >= k)
    start = (s == k).astype(float)
    Xs = sm.add_constant(d.severity.to_numpy(float))
    p = np.full(len(d), np.nan)
    p[elig] = sm.Logit(start[elig], Xs[elig]).fit(disp=0, tol=1e-12).predict(Xs[elig])
    nv = k < s
    wnv = pd.Series(np.where(k < grace, 1 / (1 - p), 1.0)[nv]).groupby(d.id[nv].to_numpy()).cumprod().to_numpy()
    gr = (s < grace - 1) | (s == grace - 1) | (k < grace - 1)
    p_last = pd.Series(p[(k == grace - 1) & (s == grace - 1)], index=d.id[(k == grace - 1) & (s == grace - 1)])
    wgr = np.where((s == grace - 1) & (k >= grace - 1), 1 / d.id.map(p_last).to_numpy(), 1.0)[gr]

    def risk(mask, w):
        kk, yy = k[mask], d.event.to_numpy(float)[mask]
        h = [np.sum(w[kk == j] * yy[kk == j]) / np.sum(w[kk == j]) if np.any(kk == j) else 0.0 for j in range(intervals)]
        return 1 - np.prod(1 - np.array(h))

    return risk(gr, wgr), risk(nv, wnv)


def target_trial(measured, truths):
    name = "target-trial-clone-censor-weight"
    regenerates(name)
    m = load(name)
    truth = m.truth()
    truths[name] = {"risk_difference": truth["risk_difference"]}
    # An independent Monte Carlo run of the truth (different seed): agreement within its sampling error.
    m.TRUTH_SEED = 12345
    assert abs(m.truth()["risk_difference"] - truth["risk_difference"]) < 0.006
    committed = run_r_many(name, [list(m.rows())])[0]
    rg, rn = ccw_numpy(pd.DataFrame(list(m.rows())))
    assert abs(rg - committed["risk_grace"]) < 1e-8 and abs(rn - committed["risk_never"]) < 1e-8, \
        "numpy clone-censor-weight disagrees with r.R"
    runs = run_r_many(name, [list(m.rows(s)) for s in SEEDS])
    t = truth["risk_difference"]
    assert np.mean([r["risk_difference_unweighted"] for r in runs]) - t > 0.02, "unweighted should be confounded"
    assert np.mean([r["risk_difference_ever_vs_never"] for r in runs]) - t < -0.01, "immortal time should favour treatment"
    measured[name] = {"risk_difference": miss99([r["risk_difference"] for r in runs], t),
                      "mean_bias_weighted": float(np.mean([r["risk_difference"] for r in runs]) - t),
                      "mean_bias_unweighted": float(np.mean([r["risk_difference_unweighted"] for r in runs]) - t),
                      "mean_bias_ever_vs_never": float(np.mean([r["risk_difference_ever_vs_never"] for r in runs]) - t)}


def check_tolerances(measured, truths):
    for name, truth in truths.items():
        exp = json.loads((ROOT / "lib" / name / "expected.json").read_text())
        for key, value in truth.items():
            assert abs(exp["truth"][key] - value) < 1e-12, f"{name} {key}: expected.json {exp['truth'][key]} != {value}"
            tol = float(exp["recovery"].get("tolerance_by_key", {}).get(key, exp["recovery"]["tolerance_estimate"]))
            assert measured[name][key] < tol, f"{name} {key}: 99th percentile miss {measured[name][key]:.4f} >= {tol}"


if __name__ == "__main__":
    measured, truths = {}, {}
    sccs(measured, truths)
    ncc(measured, truths)
    target_trial(measured, truths)
    if "--measure" in sys.argv:
        print(json.dumps({"measured": measured, "truths": truths}, indent=2))
        sys.exit(0)
    check_tolerances(measured, truths)
    print(json.dumps({"entries": len(measured), "seeds": len(SEEDS), "result": "pass"}, indent=2))
