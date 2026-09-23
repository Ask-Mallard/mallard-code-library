"""Controls and fixed-seed calibration for batch D-2.

fine-gray-subdistribution (R only), andersen-gill-recurrent-events (R and Python),
case-cohort-weighted-cox (R only). The two R-only entries are checked by R: a second R implementation
as the oracle, and recovery measured by fitting 100 generated fixtures in one Rscript call. Every entry
also gets a negative control and generator-derived truths. Run with --measure to print the numbers.

Locally, set MALLARD_R_LIBS to a library directory to put it first on R's library path; CI installs the
packages into R's default library and needs nothing.
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
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check import parse_harness_block  # noqa: E402

warnings.filterwarnings("ignore")
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


def r_prelude():
    libs = os.environ.get("MALLARD_R_LIBS")
    return f'.libPaths(c("{libs}", .libPaths())); ' if libs else ""


def run_r(code, cwd=None):
    return subprocess.check_output(["Rscript", "-e", r_prelude() + code], cwd=cwd, text=True,
                                   stderr=subprocess.DEVNULL)


def r_output(name):
    return parse_harness_block(run_r('source("r.R")', cwd=ROOT / "lib" / name))


def python_output(name):
    return parse_harness_block(subprocess.check_output(
        [sys.executable, "python.py"], cwd=ROOT / "lib" / name, text=True, stderr=subprocess.DEVNULL))


def miss99(estimates, truth):
    return rank_quantile([abs(e - truth) for e in estimates])


def fit_seeds_in_r(module, fit_code):
    """Write seeds 0..99 as CSVs and fit them all in one Rscript call; fit_code maps `d` to a number."""
    with tempfile.TemporaryDirectory() as tmp:
        for s in SEEDS:
            write_csv(Path(tmp) / f"seed_{s}.csv", list(module.rows(s)))
        out = run_r(f'for (s in 0:{len(SEEDS) - 1}) {{ d <- read.csv(file.path("{tmp}", paste0("seed_", s, ".csv"))); '
                    f'cat(sprintf("%.12f\\n", {{ {fit_code} }})) }}')
    return [float(v) for v in out.split()]


def cox_1d(entry, exit_, event, x):
    event_times = exit_[event == 1]
    xs = x[event == 1]
    risk = (entry[None, :] < event_times[:, None]) & (exit_[None, :] >= event_times[:, None])
    def nll(b):
        return -(b * xs - np.log((risk * np.exp(b * x)[None, :]).sum(axis=1))).sum()
    b = optimize.minimize_scalar(nll, bounds=(-5, 5), method="bounded", options={"xatol": 1e-12}).x
    w = risk * np.exp(b * x)[None, :]
    mean = (w * x).sum(1) / w.sum(1)
    info = ((w * x ** 2).sum(1) / w.sum(1) - mean ** 2).sum()
    return b, info, risk, mean


# --- Fine-Gray (R only) --------------------------------------------------------------------------

FG_FIT = ('suppressPackageStartupMessages(library(cmprsk)); '
          'crr(d$time, d$status, cov1 = cbind(treated = d$treated), failcode = 1, cencode = 0)$coef[1]')


def fine_gray(measured):
    name = "fine-gray-subdistribution"
    regenerates(name)
    out = r_output(name)
    module = load(name)
    # Oracle: a second R implementation, finegray() + weighted coxph with the robust variance.
    other = run_r('suppressPackageStartupMessages(library(survival)); d <- read.csv("fixture.csv"); '
                  'd$id <- seq_len(nrow(d)); fg <- finegray(Surv(time, factor(status, levels = 0:2)) ~ ., '
                  'data = d[, c("id", "time", "status", "treated")], etype = "1"); '
                  'f <- coxph(Surv(fgstart, fgstop, fgstatus) ~ treated, weights = fgwt, cluster = id, data = fg); '
                  'cat(sprintf("%.12f %.12f\\n", coef(f)[1], sqrt(vcov(f)[1, 1])))', cwd=ROOT / "lib" / name).split()
    assert abs(float(other[0]) - out["log_subdistribution_hr"]) < 1e-4
    assert abs(float(other[1]) - out["log_subdistribution_hr_se"]) < 1e-3
    # Negative control: the cause-specific hazard ratio is a different quantity.
    d = pd.read_csv(ROOT / "lib" / name / "fixture.csv")
    cs_b, *_ = cox_1d(np.zeros(len(d)), d.time.to_numpy(), (d.status == 1).to_numpy().astype(int), d.treated.to_numpy(float))
    assert abs(cs_b - out["log_subdistribution_hr"]) > 0.05
    est = fit_seeds_in_r(module, FG_FIT)
    measured[name] = {"log_subdistribution_hr": miss99(est, module.LOG_SHR),
                      "finegray_minus_crr": float(other[0]) - out["log_subdistribution_hr"],
                      "cause_specific_log_hr_committed": float(cs_b)}


# --- Andersen-Gill -------------------------------------------------------------------------------

def ag_fit(rows):
    d = pd.DataFrame(rows)
    s, t, e, x = d.start.to_numpy(), d.stop.to_numpy(), d.event.to_numpy(), d.treated.to_numpy(float)
    b, info, risk, mean = cox_1d(s, t, e, x)
    r = np.exp(b * x)
    s0 = (risk * r).sum(1)
    own = np.zeros(len(d)); own[e == 1] = x[e == 1] - mean
    share = (risk * r[None, :] * (x[None, :] - mean[:, None]) / s0[:, None]).sum(0)
    per_person = pd.Series(own - share).groupby(d.id.to_numpy()).sum().to_numpy()
    return b, 1 / math.sqrt(info), math.sqrt((per_person ** 2).sum()) / info


def andersen_gill(measured):
    name = "andersen-gill-recurrent-events"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    rows = list(module.rows())
    b, naive, robust = ag_fit(rows)
    assert abs(out["log_rate_ratio"] - b) < 1e-6 and abs(out["log_rate_ratio_naive_se"] - naive) < 1e-6
    assert abs(out["log_rate_ratio_robust_se"] - robust) < 1e-6
    # Negative controls: the naive SE is too small, and lifelines' own robust SE with delayed entry is
    # a different (smaller) number. A standing record: if lifelines changes, this says so.
    assert naive < 0.9 * robust
    from lifelines import CoxPHFitter
    lf = CoxPHFitter().fit(pd.DataFrame(rows)[["id", "treated", "start", "stop", "event"]], duration_col="stop",
                           event_col="event", entry_col="start", cluster_col="id", robust=True)
    assert abs(lf.standard_errors_["treated"] - robust) > 0.01
    est = [ag_fit(list(module.rows(s)))[0] for s in SEEDS]
    measured[name] = {"log_rate_ratio": miss99(est, module.LOG_RATE_RATIO),
                      "naive_se_committed": naive, "robust_se_committed": robust,
                      "lifelines_robust_se_committed": float(lf.standard_errors_["treated"])}


# --- case-cohort (R only) ------------------------------------------------------------------------

CC_FIT = ('suppressPackageStartupMessages(library(survival)); '
          'cch(Surv(time, event) ~ exposed, data = d, subcoh = ~subcohort, id = ~id, cohort.size = 4000, method = "LinYing")$coefficients[1]')


def case_cohort(measured):
    name = "case-cohort-weighted-cox"
    regenerates(name)
    out = r_output(name)
    module = load(name)
    # Oracle: the full cohort (which a real study never has) gives an estimate the case-cohort one
    # should be close to, within its own sampling error.
    full = pd.DataFrame(list(module.cohort()))
    full_b, *_ = cox_1d(np.zeros(len(full)), full.time.to_numpy(), full.event.to_numpy(), full.exposed.to_numpy(float))
    assert abs(out["log_hazard_ratio"] - full_b) < 2 * out["log_hazard_ratio_se"]
    # Negative control: an unweighted Cox model on the case-cohort sample, treated as if it were the cohort.
    d = pd.read_csv(ROOT / "lib" / name / "fixture.csv")
    naive_b, naive_info, *_ = cox_1d(np.zeros(len(d)), d.time.to_numpy(), d.event.to_numpy(), d.exposed.to_numpy(float))
    assert abs(naive_b - out["log_hazard_ratio"]) > 0.1
    assert 1 / math.sqrt(naive_info) < 0.7 * out["log_hazard_ratio_se"]
    est = fit_seeds_in_r(module, CC_FIT)
    measured[name] = {"log_hazard_ratio": miss99(est, module.LOG_HR),
                      "full_cohort_log_hr_committed": float(full_b),
                      "unweighted_sample_log_hr_committed": float(naive_b),
                      "unweighted_sample_se_committed": 1 / math.sqrt(naive_info)}


def generator_truth():
    return {
        "fine-gray-subdistribution": {"log_subdistribution_hr": load("fine-gray-subdistribution").LOG_SHR},
        "andersen-gill-recurrent-events": {"log_rate_ratio": load("andersen-gill-recurrent-events").LOG_RATE_RATIO},
        "case-cohort-weighted-cox": {"log_hazard_ratio": load("case-cohort-weighted-cox").LOG_HR},
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
    fine_gray(measured)
    andersen_gill(measured)
    case_cohort(measured)
    if "--measure" in sys.argv:
        print(json.dumps(measured, indent=2))
        sys.exit(0)
    check_tolerances(measured)
    print(json.dumps({"entries": len(measured), "seeds": len(SEEDS), "result": "pass"}, indent=2))
