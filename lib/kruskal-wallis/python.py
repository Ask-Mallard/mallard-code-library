# Kruskal-Wallis test across three groups, reported with each group's median
#
# The Python equivalent of r.R, using scipy's kruskal, independent of R's kruskal.test. Both apply
# the tie correction and the chi-square approximation.

import numpy as np
import pandas as pd
from scipy import stats

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any()

# PINNED: groups split by arm label, matching r.R's factor(arm).
groups = [d.y[d.arm == a].to_numpy() for a in (1, 2, 3)]
fit = stats.kruskal(*groups)
medians = [float(np.median(g)) for g in groups]

print(f"H = {fit.statistic:.3f} on {len(groups) - 1} df, p = {fit.pvalue:.4g}")
print("medians", [round(m, 3) for m in medians])

print("\n--- HARNESS ---")
print(f"h_statistic={fit.statistic:.10f}")
print(f"df={len(groups) - 1}")
print(f"p_value={fit.pvalue:.10e}")
for a, m in zip((1, 2, 3), medians):
    print(f"median_arm{a}={m:.10f}")
print(f"n={len(d)}")
