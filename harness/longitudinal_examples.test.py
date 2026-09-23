"""Controls and fixed-seed calibration for batch E-1.

linear-mixed-model-repeated-measures, cluster-level-summary-analysis, crossover-trial-two-period. Oracles
are closed forms: in a balanced random-intercept design the REML variance components are the ANOVA
estimators, and the REML estimate and standard error of a within-patient effect equal the within
(fixed-effects) regression's. Each entry also gets a negative control, generator-derived truths, and
recovery over 100 seeds. Run with --measure to print the numbers.
"""
import csv
import importlib.util
import io
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check import parse_harness_block  # noqa: E402

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


# --- linear mixed model --------------------------------------------------------------------------

def lmm_closed_form(d):
    """Within regression (subject intercepts) for the within-patient effects, plus ANOVA/REML components."""
    n, m = d.id.nunique(), d.groupby("id").size().iloc[0]
    ids = pd.get_dummies(d.id).to_numpy(float)
    X = np.column_stack([ids, d.visit, d.treated * d.visit])
    beta, *_ = np.linalg.lstsq(X, d.y.to_numpy(), rcond=None)
    resid = d.y.to_numpy() - X @ beta
    s2e = resid @ resid / (len(d) - X.shape[1])
    cov = s2e * np.linalg.inv(X.T @ X)
    g = d.groupby("id").agg(y=("y", "mean"), t=("treated", "first"))
    Xb = np.column_stack([np.ones(n), g.t])
    rb = g.y.to_numpy() - Xb @ np.linalg.lstsq(Xb, g.y.to_numpy(), rcond=None)[0]
    s2u = rb @ rb / (n - 2) - s2e / m
    return {"interaction": beta[-1], "interaction_se": math.sqrt(cov[-1, -1]), "time_slope": beta[-2],
            "sd_patient": math.sqrt(s2u), "sd_residual": math.sqrt(s2e)}


def mixed_model(measured):
    name = "linear-mixed-model-repeated-measures"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    d = pd.DataFrame(list(module.rows()))
    cf = lmm_closed_form(d)
    for key in ("interaction", "interaction_se", "time_slope"):
        assert abs(out[key] - cf[key]) < 1e-6, f"{key}: {out[key]} vs closed form {cf[key]}"
    for key in ("sd_patient", "sd_residual"):
        assert abs(out[key] - cf[key]) < 1e-6, f"{key}: {out[key]} vs closed form {cf[key]}"
    # Negative control: ordinary least squares on the 1200 rows ignores the patient structure and gives
    # a different standard error for the interaction.
    X = np.column_stack([np.ones(len(d)), d.treated, d.visit, d.treated * d.visit])
    beta, *_ = np.linalg.lstsq(X, d.y.to_numpy(), rcond=None)
    r = d.y.to_numpy() - X @ beta
    ols_se = math.sqrt(r @ r / (len(d) - 4) * np.linalg.inv(X.T @ X)[3, 3])
    assert abs(ols_se - out["interaction_se"]) > 0.02
    per_seed = [lmm_closed_form(pd.DataFrame(list(module.rows(s)))) for s in SEEDS]
    truth = generator_truth()[name]
    measured[name] = {k: miss99([p[k] for p in per_seed], truth[k]) for k in truth}
    measured[name]["ols_interaction_se_committed"] = ols_se


# --- cluster-level summary -----------------------------------------------------------------------

def cluster_summary(d):
    c = d.groupby(["cluster", "arm"], as_index=False).infection.mean()
    return c.infection[c.arm == 1].mean() - c.infection[c.arm == 0].mean(), c


def cluster_level(measured):
    name = "cluster-level-summary-analysis"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    d = pd.DataFrame(list(module.rows()))
    diff, c = cluster_summary(d)
    a, b = c.infection[c.arm == 1].to_numpy(), c.infection[c.arm == 0].to_numpy()
    sp = math.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    half = stats.t.ppf(0.975, len(a) + len(b) - 2) * sp * math.sqrt(1 / len(a) + 1 / len(b))
    assert abs(out["difference"] - diff) < 1e-9 and abs(out["difference_lcl"] - (diff - half)) < 1e-9
    # Negative control: patients analysed as independent give a much narrower interval.
    p1, p0 = d.infection[d.arm == 1].mean(), d.infection[d.arm == 0].mean()
    n1, n0 = (d.arm == 1).sum(), (d.arm == 0).sum()
    naive_half = 1.96 * math.sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
    assert naive_half < 0.8 * half
    est = [cluster_summary(pd.DataFrame(list(module.rows(s))))[0] for s in SEEDS]
    measured[name] = {"difference": miss99(est, module.TRUE_DIFFERENCE),
                      "cluster_half_width_committed": half, "naive_half_width_committed": naive_half}


# --- crossover -----------------------------------------------------------------------------------

def crossover_estimates(d):
    w = d.pivot(index=["id", "sequence"], columns="period", values="y").reset_index()
    w["diff"] = w[1] - w[2]
    adjusted = (w["diff"][w.sequence == "BA"].mean() - w["diff"][w.sequence == "AB"].mean()) / 2
    b = d[d.treatment == "B"].set_index("id").y
    a = d[d.treatment == "A"].set_index("id").y
    naive = (b - a).mean()
    return adjusted, naive, w


def crossover(measured):
    name = "crossover-trial-two-period"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    d = pd.DataFrame(list(module.rows()))
    adjusted, naive, w = crossover_estimates(d)
    # Oracle: the same effect from a regression of the period difference on a sequence indicator.
    X = np.column_stack([np.ones(len(w)), (w.sequence == "BA").astype(float)])
    coef, *_ = np.linalg.lstsq(X, w["diff"].to_numpy(), rcond=None)
    assert abs(out["treatment_effect"] - coef[1] / 2) < 1e-9 and abs(out["treatment_effect"] - adjusted) < 1e-9
    # Negative control: the paired comparison ignoring period is biased by the sequence imbalance.
    assert abs(naive - adjusted) > 0.3
    per_seed = [crossover_estimates(pd.DataFrame(list(module.rows(s)))) for s in SEEDS]
    measured[name] = {"treatment_effect": miss99([p[0] for p in per_seed], module.TREATMENT_EFFECT),
                      "naive_paired_committed": float(naive),
                      "naive_paired_bias_mean": float(np.mean([p[1] - module.TREATMENT_EFFECT for p in per_seed]))}


def generator_truth():
    lm = load("linear-mixed-model-repeated-measures")
    return {
        "linear-mixed-model-repeated-measures": {"interaction": lm.INTERACTION, "time_slope": lm.TIME_SLOPE,
                                                 "sd_patient": lm.SD_PATIENT, "sd_residual": lm.SD_NOISE},
        "cluster-level-summary-analysis": {"difference": load("cluster-level-summary-analysis").TRUE_DIFFERENCE},
        "crossover-trial-two-period": {"treatment_effect": load("crossover-trial-two-period").TREATMENT_EFFECT},
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
    mixed_model(measured)
    cluster_level(measured)
    crossover(measured)
    if "--measure" in sys.argv:
        print(json.dumps(measured, indent=2))
        sys.exit(0)
    check_tolerances(measured)
    print(json.dumps({"entries": len(measured), "seeds": len(SEEDS), "result": "pass"}, indent=2))
