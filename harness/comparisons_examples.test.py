"""Controls and fixed-seed calibration for the unadjusted comparison entries (batch B-1).

paired-t-test, one-way-anova-tukey, mann-whitney-hodges-lehmann,
wilcoxon-signed-rank-hodges-lehmann, kruskal-wallis. For each, this file regenerates the committed
fixture byte for byte, checks the Python output against an oracle written here independently, runs a
negative control for the trap the entry pins, and measures the recovery tolerance over 100 seeds.
Run with --measure to print the calibration numbers; without it, the measured 99th percentile misses
must sit inside the tolerances recorded in each expected.json.
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
        [sys.executable, "python.py"], cwd=ROOT / "lib" / name, text=True))


def miss99(estimates, truth):
    return float(np.quantile([abs(e - truth) for e in estimates], 0.99))


# --- paired t --------------------------------------------------------------------------------

def paired(measured):
    name = "paired-t-test"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    rows = list(module.rows())
    diff = np.array([r["after"] - r["before"] for r in rows])
    n = len(diff)
    se = diff.std(ddof=1) / math.sqrt(n)
    t = stats.t.ppf(0.975, n - 1)
    assert abs(out["mean_change_lcl"] - (diff.mean() - t * se)) < 1e-9
    assert abs(out["mean_change_ucl"] - (diff.mean() + t * se)) < 1e-9
    # Negative control: the unpaired (Welch) interval on the same columns is far wider.
    welch = stats.ttest_ind([r["after"] for r in rows], [r["before"] for r in rows], equal_var=False)
    welch_ci = welch.confidence_interval(0.95)
    ratio = (welch_ci.high - welch_ci.low) / (out["mean_change_ucl"] - out["mean_change_lcl"])
    assert ratio > 1.3, f"pairing barely matters on this fixture (width ratio {ratio:.2f})"
    per_seed = [np.array([r["after"] - r["before"] for r in module.rows(s)]) for s in SEEDS]
    measured[name] = {"mean_change": miss99([x.mean() for x in per_seed], module.MEAN_CHANGE),
                      "sd_change": miss99([x.std(ddof=1) for x in per_seed], module.SD_CHANGE),
                      "welch_width_ratio": ratio}


# --- one-way ANOVA with Tukey -----------------------------------------------------------------

def anova_parts(rows):
    groups = {a: np.array([r["y"] for r in rows if r["arm"] == a]) for a in (1, 2, 3)}
    n, k = sum(len(g) for g in groups.values()), 3
    mse = sum(((g - g.mean()) ** 2).sum() for g in groups.values()) / (n - k)
    return groups, n, k, mse


def anova(measured):
    name = "one-way-anova-tukey"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    rows = list(module.rows())
    groups, n, k, mse = anova_parts(rows)
    # Oracle: Tukey bounds from the studentized range quantile and the pooled variance.
    q = stats.studentized_range.ppf(0.95, k, n - k)
    for i, j in ((2, 1), (3, 1), (3, 2)):
        diff = groups[i].mean() - groups[j].mean()
        half = q * math.sqrt(mse / 2 * (1 / len(groups[i]) + 1 / len(groups[j])))
        key = f"diff_{i}_{j}"
        assert abs(out[key] - diff) < 1e-9
        assert abs(out[f"{key}_lcl"] - (diff - half)) < 1e-6 and abs(out[f"{key}_ucl"] - (diff + half)) < 1e-6
    # Negative control: arm read as NUMERIC fits one slope on 1 df, a different model and F.
    arm = np.array([r["arm"] for r in rows], dtype=float)
    y = np.array([r["y"] for r in rows])
    slope, intercept = np.polyfit(arm, y, 1)
    rss_linear = ((y - (intercept + slope * arm)) ** 2).sum()
    rss_total = ((y - y.mean()) ** 2).sum()
    f_linear = (rss_total - rss_linear) / 1 / (rss_linear / (n - 2))
    assert abs(f_linear - out["f_statistic"]) > 1.0, "numeric and factor arm are indistinguishable"
    per_seed = []
    for s in SEEDS:
        g, _, _, m = anova_parts(list(module.rows(s)))
        per_seed.append({"diff_2_1": g[2].mean() - g[1].mean(), "diff_3_1": g[3].mean() - g[1].mean(),
                         "diff_3_2": g[3].mean() - g[2].mean(), "pooled_sd": math.sqrt(m)})
    truth = {"diff_2_1": module.MEANS[2] - module.MEANS[1], "diff_3_1": module.MEANS[3] - module.MEANS[1],
             "diff_3_2": module.MEANS[3] - module.MEANS[2], "pooled_sd": module.SD}
    measured[name] = {k: miss99([e[k] for e in per_seed], truth[k]) for k in truth}
    measured[name]["f_numeric_arm"] = float(f_linear)


# --- exact rank distributions -----------------------------------------------------------------

def u_counts(m, n):
    poly = [1]
    for i in range(1, m + 1):
        poly = poly + [0] * n
        for k in range(len(poly) - 1, n + i - 1, -1):
            poly[k] -= poly[k - n - i]
        for k in range(i, len(poly)):
            poly[k] += poly[k - i]
    return poly[: m * n + 1]


def exact_two_sided(counts, observed):
    total = sum(counts)
    lower = sum(counts[: int(observed) + 1]) / total
    upper = sum(counts[int(observed):]) / total
    return min(1.0, 2 * min(lower, upper))


def hl_shift(rows):
    t = np.array([r["y"] for r in rows if r["group"] == 1])
    c = np.array([r["y"] for r in rows if r["group"] == 0])
    return float(np.median(np.subtract.outer(t, c)))


def mann_whitney(measured):
    name = "mann-whitney-hodges-lehmann"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    rows = list(module.rows())
    m = n = module.PER_GROUP
    counts = u_counts(m, n)
    # Oracle for the null distribution python.py builds: it must total C(m+n, m) arrangements and
    # reproduce scipy's exact p-value, which scipy computes by its own route.
    assert sum(counts) == math.comb(m + n, m)
    assert abs(exact_two_sided(counts, out["u_statistic"]) - out["p_value"]) < 1e-12
    assert abs(hl_shift(rows) - out["hl_shift"]) < 1e-9
    # Negative control: the difference in MEDIANS is not the Hodges-Lehmann shift.
    t = [r["y"] for r in rows if r["group"] == 1]
    c = [r["y"] for r in rows if r["group"] == 0]
    median_diff = float(np.median(t) - np.median(c))
    assert abs(median_diff - out["hl_shift"]) > 1e-3
    est = [hl_shift(list(module.rows(s))) for s in SEEDS]
    measured[name] = {"hl_shift": miss99(est, module.SHIFT), "median_difference_committed": median_diff}


def signrank_counts(n):
    poly = [1]
    for i in range(1, n + 1):
        poly = poly + [0] * i
        for k in range(len(poly) - 1, i - 1, -1):
            poly[k] += poly[k - i]
    return poly


def pseudomedian(rows):
    d = np.array([r["after"] - r["before"] for r in rows])
    return float(np.median(np.add.outer(d, d)[np.triu_indices(len(d))] / 2))


def signed_rank(measured):
    name = "wilcoxon-signed-rank-hodges-lehmann"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    rows = list(module.rows())
    counts = signrank_counts(module.N)
    assert sum(counts) == 2 ** module.N
    assert abs(exact_two_sided(counts, out["v_statistic"]) - out["p_value"]) < 1e-12
    assert abs(pseudomedian(rows) - out["pseudomedian_change"]) < 1e-9
    # Negative control: scipy's two-sided statistic is min(W+, W-), not R's V = W+. The fixture's
    # changes are mostly POSITIVE so W+ is the larger sum and the two conventions visibly differ;
    # with mostly negative changes they coincide and this control would prove nothing.
    d = np.array([r["after"] - r["before"] for r in rows])
    assert abs(stats.wilcoxon(d, method="exact").statistic - out["v_statistic"]) > 0.5
    est = [pseudomedian(list(module.rows(s))) for s in SEEDS]
    measured[name] = {"pseudomedian_change": miss99(est, module.CENTRE),
                      "scipy_statistic_committed": float(stats.wilcoxon(d, method="exact").statistic)}


# --- Kruskal-Wallis ----------------------------------------------------------------------------

def kruskal(measured):
    name = "kruskal-wallis"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    rows = list(module.rows())
    y = np.array([r["y"] for r in rows])
    arm = np.array([r["arm"] for r in rows])
    assert len(set(y)) == len(y), "the oracle below is the no-ties formula"
    ranks = stats.rankdata(y)
    big_n = len(y)
    h = 12 / (big_n * (big_n + 1)) * sum(ranks[arm == a].sum() ** 2 / (arm == a).sum() for a in (1, 2, 3)) \
        - 3 * (big_n + 1)
    assert abs(h - out["h_statistic"]) < 1e-9
    truth = {f"median_arm{a}": math.exp(mu) for a, mu in module.LOG_MEANS.items()}
    per_seed = []
    for s in SEEDS:
        rs = list(module.rows(s))
        per_seed.append({f"median_arm{a}": float(np.median([r["y"] for r in rs if r["arm"] == a])) for a in (1, 2, 3)})
    measured[name] = {k: miss99([e[k] for e in per_seed], truth[k]) for k in truth}


def check_tolerances(measured):
    for name, misses in measured.items():
        want = json.loads((ROOT / "lib" / name / "expected.json").read_text())
        rec = want["recovery"]
        for key, miss in misses.items():
            if key not in want["truth"]:
                continue
            tol = float(rec.get("tolerance_by_key", {}).get(key, rec["tolerance_estimate"]))
            assert miss < tol, f"{name} {key}: 99th percentile miss {miss:.4f} >= tolerance {tol}"


if __name__ == "__main__":
    measured = {}
    paired(measured)
    anova(measured)
    mann_whitney(measured)
    signed_rank(measured)
    kruskal(measured)
    if "--measure" in sys.argv:
        print(json.dumps(measured, indent=2))
        sys.exit(0)
    check_tolerances(measured)
    print(json.dumps({"entries": len(measured), "seeds": len(SEEDS), "result": "pass"}, indent=2))
