"""Mann-Whitney test with the Hodges-Lehmann shift and its exact 95% interval.

scipy gives the U statistic and the exact p-value but no Hodges-Lehmann interval, so the interval is
built here from the exact null distribution of U: the k-th smallest and k-th largest pairwise
differences, with k from that distribution's 2.5% quantile. This is the same construction R's
wilcox.test uses, written independently.
"""

import numpy as np
import pandas as pd
from scipy import stats

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any() and not d.y.duplicated().any()

treated = d.y[d.group == 1].to_numpy()
control = d.y[d.group == 0].to_numpy()
m, n = len(treated), len(control)


def u_null_counts(m, n):
    """Number of arrangements giving each U = 0..m*n under the null, from the generating function
    prod_{i=1..m} (1 - q^(n+i)) / (1 - q^i). Exact integers, no floating point."""
    poly = [1]
    for i in range(1, m + 1):
        poly = poly + [0] * n
        for k in range(len(poly) - 1, n + i - 1, -1):   # multiply by (1 - q^(n+i))
            poly[k] -= poly[k - n - i]
        for k in range(i, len(poly)):                   # divide by (1 - q^i)
            poly[k] += poly[k - i]
    return poly[: m * n + 1]


counts = u_null_counts(m, n)
total = sum(counts)


def q_u(p):
    """Smallest q with P(U <= q) >= p, R's qwilcox."""
    running = 0
    for q, c in enumerate(counts):
        running += c
        if running * 1.0 / total >= p:
            return q
    return m * n


# PINNED: method="exact" and treated FIRST. mannwhitneyu's default method is "auto", which switches
# to a normal approximation above 8 per group; the sign of the shift follows argument order.
test = stats.mannwhitneyu(treated, control, method="exact", alternative="two-sided")
diffs = np.sort(np.subtract.outer(treated, control).ravel())
qu = max(q_u(0.025), 1)
ql = m * n - qu
lcl, ucl = diffs[qu - 1], diffs[ql]
shift = float(np.median(diffs))

print(f"U = {test.statistic:.0f}, exact two-sided p = {test.pvalue:.4g}")
print(f"Hodges-Lehmann shift {shift:.3f}, exact 95% {lcl:.3f} to {ucl:.3f}")

print("\n--- HARNESS ---")
print(f"hl_shift={shift:.10f}")
print(f"hl_shift_lcl={lcl:.10f}")
print(f"hl_shift_ucl={ucl:.10f}")
print(f"u_statistic={test.statistic:.10f}")
print(f"p_value={test.pvalue:.10e}")
print(f"n={len(d)}")
