"""Controls and fixed-seed calibration for the descriptive entries.

Four entries: single-proportion-wilson-clopper-pearson, incidence-rate-exact-poisson,
clustered-prevalence-cluster-sample, descriptive-summary-table-one. For each, this file:

  1. regenerates the committed fixture from its generator and requires it byte for byte;
  2. checks the committed Python output against an ORACLE written here independently, from the
     textbook formula, so R and Python agreeing with each other is not the only evidence;
  3. runs a NEGATIVE CONTROL for the trap the entry pins, showing the wrong method gives a number
     the checks would catch (Wald instead of Wilson, a count interval instead of a rate interval,
     the average of clinic percentages instead of the patient-level prevalence, ddof=0 instead of
     n - 1);
  4. measures the recovery tolerance over 100 seeds and requires the 99th percentile miss to sit
     inside the tolerance recorded in expected.json. Run with --measure to print the numbers.

A tolerance chosen by assertion is a check calibrated against nothing (README.md records the 1e-6
that failed a correct pair). The numbers in each expected.json came from this file.
"""
import csv
import importlib.util
import io
import json
import math
import statistics
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check import parse_harness_block  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
SEEDS = range(100)
Z = statistics.NormalDist().inv_cdf(0.975)


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
        [sys.executable, "python.py"], cwd=ROOT / "lib" / name, text=True))


def expected(name):
    return json.loads((ROOT / "lib" / name / "expected.json").read_text())


def tolerance(name, key):
    rec = expected(name)["recovery"]
    return float(rec.get("tolerance_by_key", {}).get(key, rec["tolerance_estimate"]))


def miss99(estimates, truth):
    return float(np.quantile([abs(e - truth) for e in estimates], 0.99))


# --- single proportion ----------------------------------------------------------------------

def wilson(x, n):
    p = x / n
    centre = (p + Z * Z / (2 * n)) / (1 + Z * Z / n)
    half = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / (1 + Z * Z / n)
    return centre - half, centre + half


def single_proportion(measured):
    name = "single-proportion-wilson-clopper-pearson"
    regenerates(name)
    out = python_output(name)
    x, n = int(out["events"]), int(out["n"])
    lo, hi = wilson(x, n)
    assert abs(out["wilson_lcl"] - lo) < 1e-9 and abs(out["wilson_ucl"] - hi) < 1e-9
    assert abs(out["cp_lcl"] - stats.beta.ppf(0.025, x, n - x + 1)) < 1e-9
    assert abs(out["cp_ucl"] - stats.beta.ppf(0.975, x + 1, n - x)) < 1e-9
    # Negative control: the Wald interval must differ by more than the agreement tolerance, or a
    # silent fall-back to Wald would pass.
    p = x / n
    wald_lo = p - Z * math.sqrt(p * (1 - p) / n)
    assert abs(wald_lo - lo) > 1e-3, "Wald and Wilson are indistinguishable on this fixture"
    module = load(name)
    estimates = [sum(r["outcome"] for r in module.rows(seed)) / module.N for seed in SEEDS]
    measured[name] = {"prevalence": miss99(estimates, module.PROPORTION),
                      "wald_minus_wilson_lower": wald_lo - lo}


# --- incidence rate -------------------------------------------------------------------------

def incidence_rate(measured):
    name = "incidence-rate-exact-poisson"
    regenerates(name)
    out = python_output(name)
    x, pt = int(out["events"]), out["person_years"]
    # Garwood's exact interval, from chi-square quantiles: a different route from both engines.
    assert abs(out["rate_lcl"] - stats.chi2.ppf(0.025, 2 * x) / 2 / pt) < 1e-9
    assert abs(out["rate_ucl"] - stats.chi2.ppf(0.975, 2 * x + 2) / 2 / pt) < 1e-9
    # Negative control: the interval for the COUNT (poisson.test's default T = 1) is not the rate's.
    assert stats.chi2.ppf(0.025, 2 * x) / 2 > out["rate_ucl"] * 100
    module = load(name)
    estimates = []
    for seed in SEEDS:
        rows = list(module.rows(seed))
        estimates.append(sum(r["events"] for r in rows) / sum(r["followup_years"] for r in rows))
    measured[name] = {"rate": miss99(estimates, module.RATE)}


# --- clustered prevalence -------------------------------------------------------------------

def clustered_prevalence(measured):
    name = "clustered-prevalence-cluster-sample"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    rows = list(module.rows())
    clusters = {}
    for r in rows:
        y, n = clusters.get(r["cluster"], (0, 0))
        clusters[r["cluster"]] = (y + r["outcome"], n + 1)
    # Oracle: the ratio estimator's with-replacement variance written from cluster totals,
    # v = m/(m-1) * sum (y_i - r n_i)^2 / (sum n_i)^2, algebraically the linearization in python.py
    # but computed a different way, from counts rather than weighted residuals.
    m = len(clusters)
    total_n = sum(n for _, n in clusters.values())
    ratio = sum(y for y, _ in clusters.values()) / total_n
    se = math.sqrt(m / (m - 1) * sum((y - ratio * n) ** 2 for y, n in clusters.values())) / total_n
    assert abs(out["prevalence"] - ratio) < 1e-9 and abs(out["prevalence_se"] - se) < 1e-9
    # Negative controls. (a) Treating patients as independent gives a smaller SE.
    srs_se = math.sqrt(ratio * (1 - ratio) / total_n)
    assert srs_se < 0.8 * se, "the design effect is too small to show on this fixture"
    # (b) The average of the clinic percentages estimates a different quantity and fails recovery.
    # Asserted against the recorded tolerance in check_tolerances, once expected.json exists.
    cluster_mean = statistics.mean(y / n for y, n in clusters.values())
    estimates, naive = [], []
    for seed in SEEDS:
        cl = {}
        for r in module.rows(seed):
            y, n = cl.get(r["cluster"], (0, 0))
            cl[r["cluster"]] = (y + r["outcome"], n + 1)
        estimates.append(sum(y for y, _ in cl.values()) / sum(n for _, n in cl.values()))
        naive.append(statistics.mean(y / n for y, n in cl.values()))
    measured[name] = {"prevalence": miss99(estimates, module.TRUE_PREVALENCE),
                      "design_effect_se_ratio": se / srs_se,
                      "cluster_mean_miss_committed": abs(cluster_mean - module.TRUE_PREVALENCE),
                      "cluster_mean_miss_median": float(np.median(
                          [abs(v - module.TRUE_PREVALENCE) for v in naive]))}


# --- Table 1 --------------------------------------------------------------------------------

def table_one_truth(module):
    z25 = statistics.NormalDist().inv_cdf(0.25)
    p1, p0 = module.DIABETES[1], module.DIABETES[0]
    spread = (module.AGE_MEAN[1] - module.AGE_MEAN[0]) / 2
    return {
        "age_mean": (module.AGE_MEAN[0] + module.AGE_MEAN[1]) / 2,
        "age_sd": math.sqrt(module.AGE_SD ** 2 + spread ** 2),
        "age_mean_exposed": module.AGE_MEAN[1], "age_mean_unexposed": module.AGE_MEAN[0],
        "age_sd_exposed": module.AGE_SD, "age_sd_unexposed": module.AGE_SD,
        "age_smd": (module.AGE_MEAN[1] - module.AGE_MEAN[0]) / module.AGE_SD,
        "los_q1": math.exp(module.LOS_MU + z25 * module.LOS_SIGMA),
        "los_median": math.exp(module.LOS_MU),
        "los_q3": math.exp(module.LOS_MU - z25 * module.LOS_SIGMA),
        "diabetes_prop": (p0 + p1) / 2,
        "diabetes_prop_exposed": p1, "diabetes_prop_unexposed": p0,
        "diabetes_smd": (p1 - p0) / math.sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2),
        **{f"smoking_{level}_prop": p for level, p in module.SMOKING},
    }


def table_one_estimates(rows):
    e1 = [r for r in rows if r["exposed"] == 1]
    e0 = [r for r in rows if r["exposed"] == 0]

    def col(rs, k):
        return [r[k] for r in rs]

    def smd(a, b):
        return (statistics.mean(a) - statistics.mean(b)) / math.sqrt(
            (statistics.stdev(a) ** 2 + statistics.stdev(b) ** 2) / 2)

    def smd_bin(a, b):
        p1, p0 = statistics.mean(a), statistics.mean(b)
        return (p1 - p0) / math.sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)

    # statistics.quantiles(method="inclusive") is R's type 7, an implementation separate from both
    # engines. statistics.stdev divides by n - 1.
    q1, med, q3 = statistics.quantiles(col(rows, "los_days"), n=4, method="inclusive")
    out = {
        "age_mean": statistics.mean(col(rows, "age")), "age_sd": statistics.stdev(col(rows, "age")),
        "age_mean_exposed": statistics.mean(col(e1, "age")),
        "age_mean_unexposed": statistics.mean(col(e0, "age")),
        "age_sd_exposed": statistics.stdev(col(e1, "age")),
        "age_sd_unexposed": statistics.stdev(col(e0, "age")),
        "age_smd": smd(col(e1, "age"), col(e0, "age")),
        "los_q1": q1, "los_median": med, "los_q3": q3,
        "diabetes_prop": statistics.mean(col(rows, "diabetes")),
        "diabetes_prop_exposed": statistics.mean(col(e1, "diabetes")),
        "diabetes_prop_unexposed": statistics.mean(col(e0, "diabetes")),
        "diabetes_smd": smd_bin(col(e1, "diabetes"), col(e0, "diabetes")),
    }
    for level in ("never", "former", "current"):
        out[f"smoking_{level}_prop"] = sum(r["smoking"] == level for r in rows) / len(rows)
    return out


def table_one(measured):
    name = "descriptive-summary-table-one"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    rows = list(module.rows())
    oracle = table_one_estimates(rows)
    for key, value in oracle.items():
        assert abs(out[key] - value) < 1e-9, f"{key}: python.py {out[key]} vs oracle {value}"
    # Negative controls: the n divisor and a different quartile definition must each move a number
    # by more than the agreement tolerance, or dropping the pin would go unnoticed.
    ages = [r["age"] for r in rows]
    assert abs(statistics.pstdev(ages) - out["age_sd"]) > 1e-3, "ddof=0 is indistinguishable"
    los = [r["los_days"] for r in rows]
    exclusive = statistics.quantiles(los, n=4, method="exclusive")
    assert max(abs(exclusive[0] - out["los_q1"]), abs(exclusive[2] - out["los_q3"])) > 1e-3, \
        "type 6 and type 7 quartiles are indistinguishable"
    truth = table_one_truth(module)
    per_seed = [table_one_estimates(list(module.rows(seed))) for seed in SEEDS]
    measured[name] = {k: miss99([e[k] for e in per_seed], truth[k]) for k in truth}


def check_tolerances(measured):
    # The wrong estimand must fail recovery on the committed fixture, or the check cannot see it.
    cp = measured["clustered-prevalence-cluster-sample"]
    assert cp["cluster_mean_miss_committed"] > tolerance("clustered-prevalence-cluster-sample", "prevalence"), \
        "the average of clinic percentages would pass recovery on this fixture"
    for name, misses in measured.items():
        want = expected(name)
        for key, miss in misses.items():
            if key not in want["truth"]:
                continue
            tol = tolerance(name, key)
            assert miss < tol, f"{name} {key}: 99th percentile miss {miss:.4f} >= tolerance {tol}"
            assert abs(want["truth"][key] - (table_one_truth(load(name))[key]
                       if name == "descriptive-summary-table-one" else want["truth"][key])) < 1e-12


if __name__ == "__main__":
    measured = {}
    single_proportion(measured)
    incidence_rate(measured)
    clustered_prevalence(measured)
    table_one(measured)
    if "--measure" in sys.argv:
        print(json.dumps(measured, indent=2))
        sys.exit(0)
    check_tolerances(measured)
    print(json.dumps({"entries": len(measured), "seeds": len(SEEDS), "result": "pass"}, indent=2))
