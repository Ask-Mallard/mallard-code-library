"""Wilcoxon signed-rank test with the Hodges-Lehmann pseudo-median change and its exact 95% interval.

scipy gives the exact p-value but no pseudo-median or interval, and its reported statistic for a
two-sided test is min(W+, W-), not the W+ (sum of positive ranks) R reports as V. So V is computed
here, and the interval is built from the exact null distribution of the signed-rank statistic: the
k-th smallest and k-th largest Walsh averages.
"""

import numpy as np
import pandas as pd
from scipy import stats

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any()

change = (d.after - d.before).to_numpy()
n = len(change)
assert (change != 0).all(), "zero changes need a stated handling rule"


def signed_rank_null_counts(n):
    """Number of sign patterns giving each W+ = 0..n(n+1)/2, from prod_{i=1..n} (1 + q^i)."""
    poly = [1]
    for i in range(1, n + 1):
        poly = poly + [0] * i
        for k in range(len(poly) - 1, i - 1, -1):
            poly[k] += poly[k - i]
    return poly


counts = signed_rank_null_counts(n)
total = sum(counts)


def q_signrank(p):
    """Smallest q with P(W+ <= q) >= p, R's qsignrank."""
    running = 0
    for q, c in enumerate(counts):
        running += c
        if running * 1.0 / total >= p:
            return q
    return len(counts) - 1


# PINNED: V is the sum of the ranks of the POSITIVE changes, R's convention. scipy.stats.wilcoxon's
# two-sided statistic is min(W+, W-); its p-value is used, not its statistic. method="exact" pinned.
ranks = stats.rankdata(np.abs(change))
v = float(ranks[change > 0].sum())
p = stats.wilcoxon(change, method="exact", alternative="two-sided").pvalue

walsh = np.sort(np.add.outer(change, change)[np.triu_indices(n)] / 2)
qu = max(q_signrank(0.025), 1)
ql = n * (n + 1) // 2 - qu
lcl, ucl = walsh[qu - 1], walsh[ql]
pseudomedian = float(np.median(walsh))

print(f"V = {v:.0f}, exact two-sided p = {p:.4g}")
print(f"pseudo-median change {pseudomedian:.3f}, exact 95% {lcl:.3f} to {ucl:.3f}")

print("\n--- HARNESS ---")
print(f"pseudomedian_change={pseudomedian:.10f}")
print(f"pseudomedian_change_lcl={lcl:.10f}")
print(f"pseudomedian_change_ucl={ucl:.10f}")
print(f"v_statistic={v:.10f}")
print(f"p_value={p:.10e}")
print(f"n={n}")
