"""Controls and fixed-seed calibration for the categorical and correlation entries.

fisher-exact-conditional-odds-ratio, mcnemar-paired-binary, risk-difference-ratio-nnt,
correlation-pearson-spearman, cochran-armitage-trend. For each: byte-for-byte fixture regeneration,
an oracle written here independently, a negative control for the pinned trap, and recovery measured
over 100 seeds. Run with --measure to print the calibration numbers.
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
from scipy import optimize, stats

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
Z = stats.norm.ppf(0.975)


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
    return rank_quantile([abs(e - truth) for e in estimates])


def cells(rows, row_key, col_key):
    a = sum(1 for r in rows if r[row_key] == 1 and r[col_key] == 1)
    b = sum(1 for r in rows if r[row_key] == 1 and r[col_key] == 0)
    c = sum(1 for r in rows if r[row_key] == 0 and r[col_key] == 1)
    d = sum(1 for r in rows if r[row_key] == 0 and r[col_key] == 0)
    return a, b, c, d


# --- Fisher --------------------------------------------------------------------------------------

def conditional_log_or(a, b, c, d):
    """Conditional MLE: the odds ratio at which the Fisher noncentral hypergeometric mean of the
    top-left cell equals its observed count. A different code path from contingency.odds_ratio."""
    if 0 in (a, b, c, d):
        return math.inf
    total, row1, col1 = a + b + c + d, a + b, a + c
    f = lambda lo: stats.nchypergeom_fisher.mean(total, row1, col1, math.exp(lo)) - a
    return optimize.brentq(f, -20, 20, xtol=1e-12)


def fisher(measured):
    name = "fisher-exact-conditional-odds-ratio"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    a, b, c, d = cells(list(module.rows()), "exposed", "outcome")
    assert abs(conditional_log_or(a, b, c, d) - out["log_odds_ratio"]) < 1e-6
    # Negative control: the sample odds ratio (fisher_exact's statistic) is not the conditional MLE.
    # The two converge as tables grow (0.044 apart here), so the bar is ten times the agreement
    # tolerance: enough that reporting fisher_exact's statistic would fail agreement.
    assert abs(out["sample_odds_ratio"] - out["odds_ratio"]) > 1e-3
    est = [conditional_log_or(*cells(list(module.rows(s)), "exposed", "outcome")) for s in SEEDS]
    infinite = sum(1 for e in est if math.isinf(e))
    measured[name] = {"log_odds_ratio": miss99([e for e in est if not math.isinf(e)], module.LOG_OR),
                      "seeds_with_a_zero_cell": infinite}


# --- McNemar -------------------------------------------------------------------------------------

def mcnemar_check(measured):
    name = "mcnemar-paired-binary"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    rows = list(module.rows())
    b = sum(1 for r in rows if r["before"] == 1 and r["after"] == 0)
    c = sum(1 for r in rows if r["before"] == 0 and r["after"] == 1)
    assert abs(out["mcnemar_chi2"] - (b - c) ** 2 / (b + c)) < 1e-9
    assert abs(out["exact_p_value"] - min(1.0, 2 * stats.binom.cdf(min(b, c), b + c, 0.5))) < 1e-12
    # Negative controls: the continuity-corrected statistic, and an UNPAIRED chi-square on the same
    # columns, are both different numbers.
    corrected = (abs(b - c) - 1) ** 2 / (b + c)
    assert abs(corrected - out["mcnemar_chi2"]) > 0.1
    n = len(rows)
    before, after = sum(r["before"] for r in rows), sum(r["after"] for r in rows)
    unpaired = stats.chi2_contingency([[before, n - before], [after, n - after]], correction=False)[0]
    assert abs(unpaired - out["mcnemar_chi2"]) > 0.1
    est = []
    for s in SEEDS:
        rs = list(module.rows(s))
        est.append((sum(r["after"] for r in rs) - sum(r["before"] for r in rs)) / len(rs))
    measured[name] = {"change": miss99(est, module.TRUE_CHANGE), "unpaired_chi2_committed": float(unpaired)}


# --- risk difference, ratio, NNT -----------------------------------------------------------------

def wilson(x, n):
    p = x / n
    centre = (p + Z * Z / (2 * n)) / (1 + Z * Z / n)
    half = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / (1 + Z * Z / n)
    return centre - half, centre + half


def risk_effects(measured):
    name = "risk-difference-ratio-nnt"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    rows = list(module.rows())
    x1 = sum(r["outcome"] for r in rows if r["treated"] == 1); n1 = sum(r["treated"] for r in rows)
    x0 = sum(r["outcome"] for r in rows if r["treated"] == 0); n0 = len(rows) - n1
    p1, p0 = x1 / n1, x0 / n0
    l1, u1 = wilson(x1, n1)
    l0, u0 = wilson(x0, n0)
    rd = p1 - p0
    assert abs(out["risk_difference_lcl"] - (rd - math.sqrt((p1 - l1) ** 2 + (u0 - p0) ** 2))) < 1e-9
    assert abs(out["risk_difference_ucl"] - (rd + math.sqrt((u1 - p1) ** 2 + (p0 - l0) ** 2))) < 1e-9
    se = math.sqrt(1 / x1 - 1 / n1 + 1 / x0 - 1 / n0)
    assert abs(math.log(out["risk_ratio_lcl"]) - (math.log(p1 / p0) - Z * se)) < 1e-9
    # Negative control: the Wald RD interval differs from Newcombe's.
    wald_lcl = rd - Z * math.sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
    assert abs(wald_lcl - out["risk_difference_lcl"]) > 1e-4
    rds, lrrs, crossing = [], [], 0
    for s in SEEDS:
        rs = list(module.rows(s))
        a1 = sum(r["outcome"] for r in rs if r["treated"] == 1)
        a0 = sum(r["outcome"] for r in rs if r["treated"] == 0)
        q1, q0 = a1 / module.PER_ARM, a0 / module.PER_ARM
        rds.append(q1 - q0)
        lrrs.append(math.log(q1 / q0))
        w1, v1 = wilson(a1, module.PER_ARM)
        w0, v0 = wilson(a0, module.PER_ARM)
        if (q1 - q0) + math.sqrt((v1 - q1) ** 2 + (q0 - w0) ** 2) >= 0:
            crossing += 1
    measured[name] = {"risk_difference": miss99(rds, module.RISK_DIFFERENCE),
                      "log_risk_ratio": miss99(lrrs, module.LOG_RISK_RATIO),
                      "wald_minus_newcombe_lcl": wald_lcl - out["risk_difference_lcl"],
                      "seeds_whose_rd_interval_crosses_zero": crossing}


# --- correlation ---------------------------------------------------------------------------------

def correlations(rows):
    x = np.array([r["bmi"] for r in rows]); y = np.array([r["sbp"] for r in rows])
    return float(np.corrcoef(x, y)[0, 1]), float(np.corrcoef(stats.rankdata(x), stats.rankdata(y))[0, 1])


def correlation(measured):
    name = "correlation-pearson-spearman"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    rows = list(module.rows())
    r, rho = correlations(rows)
    n = len(rows)
    assert abs(out["pearson_r"] - r) < 1e-9 and abs(out["spearman_rho"] - rho) < 1e-9
    assert abs(out["pearson_lcl"] - math.tanh(math.atanh(r) - Z / math.sqrt(n - 3))) < 1e-9
    # Negative control: Pearson's standard error on Spearman gives a narrower interval.
    naive = math.tanh(math.atanh(rho) - Z / math.sqrt(n - 3))
    assert abs(naive - out["spearman_lcl"]) > 1e-3
    per_seed = [correlations(list(module.rows(s))) for s in SEEDS]
    measured[name] = {"pearson_r": miss99([p[0] for p in per_seed], module.RHO),
                      "spearman_rho": miss99([p[1] for p in per_seed], module.SPEARMAN)}


# --- Cochran-Armitage ----------------------------------------------------------------------------

def slope(rows):
    doses = sorted({r["dose"] for r in rows})
    n = np.array([sum(1 for r in rows if r["dose"] == k) for k in doses], dtype=float)
    p = np.array([sum(r["outcome"] for r in rows if r["dose"] == k) for k in doses]) / n
    s = np.array(doses, dtype=float)
    sb, pb = (n * s).sum() / n.sum(), (n * p).sum() / n.sum()
    return float((n * (s - sb) * (p - pb)).sum() / (n * (s - sb) ** 2).sum())


def trend(measured):
    name = "cochran-armitage-trend"
    regenerates(name)
    out = python_output(name)
    module = load(name)
    rows = list(module.rows())
    y = np.array([r["outcome"] for r in rows], dtype=float)
    x = np.array([r["dose"] for r in rows], dtype=float)
    r = np.corrcoef(x, y)[0, 1]
    big_n = len(rows)
    # Oracle: the Cochran-Armitage statistic is N r^2 over individuals.
    assert abs(out["trend_chi2"] - big_n * r * r) < 1e-9
    # Negative control: the linear-by-linear association statistic is (N - 1) r^2.
    assert abs((big_n - 1) * r * r - out["trend_chi2"]) > 1e-4
    est = [slope(list(module.rows(s))) for s in SEEDS]
    measured[name] = {"slope": miss99(est, module.SLOPE),
                      "linear_by_linear_minus_ca": float((big_n - 1) * r * r - out["trend_chi2"])}


def generator_truth():
    """Every truth value expected.json states, recomputed from the generator's own constants. A typed
    truth can be wrong by less than the recovery tolerance and still pass recovery; this cannot."""
    return {
        "fisher-exact-conditional-odds-ratio": {"log_odds_ratio": load("fisher-exact-conditional-odds-ratio").LOG_OR},
        "mcnemar-paired-binary": {"change": load("mcnemar-paired-binary").TRUE_CHANGE},
        "risk-difference-ratio-nnt": {"risk_difference": load("risk-difference-ratio-nnt").RISK_DIFFERENCE,
                                      "log_risk_ratio": load("risk-difference-ratio-nnt").LOG_RISK_RATIO},
        "correlation-pearson-spearman": {"pearson_r": load("correlation-pearson-spearman").RHO,
                                         "spearman_rho": load("correlation-pearson-spearman").SPEARMAN},
        "cochran-armitage-trend": {"slope": load("cochran-armitage-trend").SLOPE},
    }


def check_tolerances(measured):
    for name, truth in generator_truth().items():
        stated = json.loads((ROOT / "lib" / name / "expected.json").read_text())["truth"]
        numeric = {k: v for k, v in stated.items() if isinstance(v, (int, float))}
        assert set(numeric) == set(truth), f"{name}: expected.json truth keys {sorted(numeric)} != generator {sorted(truth)}"
        for key, value in truth.items():
            assert abs(numeric[key] - value) < 1e-12, f"{name} {key}: expected.json {numeric[key]} != generator {value}"
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
    fisher(measured)
    mcnemar_check(measured)
    risk_effects(measured)
    correlation(measured)
    trend(measured)
    if "--measure" in sys.argv:
        print(json.dumps(measured, indent=2))
        sys.exit(0)
    check_tolerances(measured)
    print(json.dumps({"entries": len(measured), "seeds": len(SEEDS), "result": "pass"}, indent=2))
