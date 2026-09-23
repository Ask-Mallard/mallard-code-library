"""Controls and fixed-seed calibration for the survival entries in batch D-1.

cox-time-varying-exposure, cause-specific-cox, restricted-mean-survival-time,
weibull-accelerated-failure-time. Every oracle here is computed from scratch in numpy/scipy (partial
likelihoods, a Kaplan-Meier estimate, a Weibull likelihood), independent of survival and lifelines. Each
entry also gets a negative control, generator-derived truths, and recovery over 100 seeds. Run with
--measure to print the calibration numbers.
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
from scipy import optimize

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check import parse_harness_block  # noqa: E402

warnings.filterwarnings("ignore")
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


def cox_1d(entry, exit_, event, x):
    """Maximize a one-covariate Cox partial likelihood (counting process, no tied event times)."""
    event_times = exit_[event == 1]
    xs = x[event == 1]
    risk = (entry[None, :] < event_times[:, None]) & (exit_[None, :] >= event_times[:, None])
    def nll(b):
        return -(b * xs - np.log((risk * np.exp(b * x)[None, :]).sum(axis=1))).sum()
    res = optimize.minimize_scalar(nll, bounds=(-5, 5), method="bounded", options={"xatol": 1e-12})
    b = res.x
    w = risk * np.exp(b * x)[None, :]
    mean = (w * x).sum(1) / w.sum(1)
    info = ((w * x ** 2).sum(1) / w.sum(1) - mean ** 2).sum()
    return b, 1 / math.sqrt(info)


# --- time-varying Cox ----------------------------------------------------------------------------

def tv_arrays(rows):
    return (np.array([r["start"] for r in rows]), np.array([r["stop"] for r in rows]),
            np.array([r["event"] for r in rows]), np.array([r["exposed"] for r in rows], dtype=float))


def time_varying(measured):
    name = "cox-time-varying-exposure"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    rows = list(module.rows())
    b, se = cox_1d(*tv_arrays(rows))
    assert abs(out["log_hazard_ratio"] - b) < 1e-6 and abs(out["log_hazard_ratio_se"] - se) < 1e-6
    # Negative control: "ever exposed" at baseline, one row per person from time 0. Immortal time.
    people = {}
    for r in rows:
        p = people.setdefault(r["id"], {"stop": 0.0, "event": 0, "ever": 0})
        p["stop"] = max(p["stop"], r["stop"]); p["event"] |= r["event"]; p["ever"] |= r["exposed"]
    arr = list(people.values())
    ever_b, _ = cox_1d(np.zeros(len(arr)), np.array([p["stop"] for p in arr]),
                       np.array([p["event"] for p in arr]), np.array([p["ever"] for p in arr], dtype=float))
    est = [cox_1d(*tv_arrays(list(module.rows(s))))[0] for s in SEEDS]
    measured[name] = {"log_hazard_ratio": miss99(est, module.LOG_HR), "ever_exposed_log_hr_committed": float(ever_b)}
    measured[name]["_ever_minus_truth"] = abs(ever_b - module.LOG_HR)


# --- cause-specific Cox --------------------------------------------------------------------------

def cs_arrays(rows, composite=False):
    t = np.array([r["time"] for r in rows])
    s = np.array([r["status"] for r in rows])
    e = (s > 0).astype(int) if composite else (s == 1).astype(int)
    return np.zeros(len(t)), t, e, np.array([r["treated"] for r in rows], dtype=float)


def cause_specific(measured):
    name = "cause-specific-cox"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    rows = list(module.rows())
    b, se = cox_1d(*cs_arrays(rows))
    assert abs(out["log_cs_hazard_ratio"] - b) < 1e-6 and abs(out["log_cs_hazard_ratio_se"] - se) < 1e-6
    # Negative control: counting deaths as events (a composite) estimates something else.
    composite_b, _ = cox_1d(*cs_arrays(rows, composite=True))
    assert abs(composite_b - out["log_cs_hazard_ratio"]) > 0.05
    est = [cox_1d(*cs_arrays(list(module.rows(s))))[0] for s in SEEDS]
    measured[name] = {"log_cs_hazard_ratio": miss99(est, module.LOG_CS_HR), "composite_log_hr_committed": float(composite_b)}


# --- RMST ----------------------------------------------------------------------------------------

def km_rmst(time, event, tau):
    """Kaplan-Meier area to tau and its Greenwood-type standard error, from the raw times."""
    ut = np.unique(time[event == 1])
    ut = ut[ut <= tau]
    n = np.array([(time >= t).sum() for t in ut], dtype=float)
    d = np.array([((time == t) & (event == 1)).sum() for t in ut], dtype=float)
    s_after = np.cumprod(1 - d / n)                      # S just after each event time
    pieces = s_after * np.diff(np.append(ut, tau))       # flat at S(t_j) until the next event time
    area = (ut[0] if len(ut) else tau) + pieces.sum()
    after = np.cumsum(pieces[::-1])[::-1]                # area from t_j to tau
    return float(area), float(np.sqrt((after ** 2 * d / (n * (n - d))).sum()))


def rmst(measured):
    name = "restricted-mean-survival-time"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    rows = list(module.rows())
    arms = {a: (np.array([r["time"] for r in rows if r["arm"] == a]), np.array([r["event"] for r in rows if r["arm"] == a]))
            for a in (0, 1)}
    (r1, s1), (r0, s0) = km_rmst(*arms[1], module.TAU), km_rmst(*arms[0], module.TAU)
    assert abs(out["rmst_treated"] - r1) < 1e-9 and abs(out["rmst_control"] - r0) < 1e-9
    assert abs(out["rmst_treated_se"] - s1) < 1e-9 and abs(out["rmst_control_se"] - s0) < 1e-9
    # Negative control: lifelines' return_variance is the variance of min(T, tau), not of the estimate.
    from lifelines import KaplanMeierFitter
    from lifelines.utils import restricted_mean_survival_time
    _, v = restricted_mean_survival_time(KaplanMeierFitter().fit(*arms[1]), t=module.TAU, return_variance=True)
    assert math.sqrt(v) > 5 * out["rmst_treated_se"]
    per_seed = []
    for s in SEEDS:
        rs = list(module.rows(s))
        a1 = km_rmst(np.array([r["time"] for r in rs if r["arm"] == 1]), np.array([r["event"] for r in rs if r["arm"] == 1]), module.TAU)[0]
        a0 = km_rmst(np.array([r["time"] for r in rs if r["arm"] == 0]), np.array([r["event"] for r in rs if r["arm"] == 0]), module.TAU)[0]
        per_seed.append((a1, a0))
    measured[name] = {"rmst_treated": miss99([p[0] for p in per_seed], module.TRUE_RMST[1]),
                      "rmst_control": miss99([p[1] for p in per_seed], module.TRUE_RMST[0]),
                      "rmst_difference": miss99([p[0] - p[1] for p in per_seed], module.TRUE_DIFFERENCE),
                      "lifelines_variance_as_se_committed": math.sqrt(v)}


# --- Weibull AFT ---------------------------------------------------------------------------------

def weibull_mle(rows):
    t = np.array([r["time"] for r in rows]); e = np.array([r["event"] for r in rows]); x = np.array([r["treated"] for r in rows], float)
    def nll(p):
        b0, b1, ls = p; s = math.exp(ls); z = (np.log(t) - b0 - b1 * x) / s
        return -(e * (z - ls - np.log(t)) - np.exp(z)).sum()
    res = optimize.minimize(nll, [float(np.log(t).mean()), 0.0, 0.0], method="BFGS", options={"gtol": 1e-10})
    b0, b1, ls = res.x
    return b1, 1 / math.exp(ls)


def weibull(measured):
    name = "weibull-accelerated-failure-time"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    rows = list(module.rows())
    b1, shape = weibull_mle(rows)
    assert abs(out["log_time_ratio"] - b1) < 1e-4 and abs(out["shape"] - shape) < 1e-4
    # Negative control: the AFT coefficient is not a log hazard ratio; they differ in sign here.
    assert out["log_time_ratio"] > 0 > out["log_hazard_ratio"]
    per_seed = [weibull_mle(list(module.rows(s))) for s in SEEDS]
    measured[name] = {"log_time_ratio": miss99([p[0] for p in per_seed], module.LOG_TIME_RATIO),
                      "shape": miss99([p[1] for p in per_seed], module.SHAPE),
                      "log_hazard_ratio": miss99([-p[0] * p[1] for p in per_seed], module.LOG_HAZARD_RATIO)}


def generator_truth():
    tv = load("cox-time-varying-exposure"); cs = load("cause-specific-cox")
    rm = load("restricted-mean-survival-time"); wb = load("weibull-accelerated-failure-time")
    return {
        "cox-time-varying-exposure": {"log_hazard_ratio": tv.LOG_HR},
        "cause-specific-cox": {"log_cs_hazard_ratio": cs.LOG_CS_HR},
        "restricted-mean-survival-time": {"rmst_treated": rm.TRUE_RMST[1], "rmst_control": rm.TRUE_RMST[0],
                                          "rmst_difference": rm.TRUE_DIFFERENCE},
        "weibull-accelerated-failure-time": {"log_time_ratio": wb.LOG_TIME_RATIO, "shape": wb.SHAPE,
                                             "log_hazard_ratio": wb.LOG_HAZARD_RATIO},
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
    # The immortal-time negative control must fail recovery by a wide margin.
    tv = measured["cox-time-varying-exposure"]
    rec = json.loads((ROOT / "lib/cox-time-varying-exposure/expected.json").read_text())["recovery"]
    assert tv["_ever_minus_truth"] > float(rec["tolerance_estimate"]), "ever-exposed coding would pass recovery"


if __name__ == "__main__":
    measured = {}
    time_varying(measured)
    cause_specific(measured)
    rmst(measured)
    weibull(measured)
    if "--measure" in sys.argv:
        print(json.dumps(measured, indent=2))
        sys.exit(0)
    check_tolerances(measured)
    print(json.dumps({"entries": len(measured), "seeds": len(SEEDS), "result": "pass"}, indent=2))
