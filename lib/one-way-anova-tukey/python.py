# One-way ANOVA across three groups, then all pairwise differences with Tukey-adjusted 95% intervals
#
# The Python equivalent of r.R: scipy's f_oneway and tukey_hsd, implementations independent of R's
# aov and TukeyHSD (scipy computes the studentized range distribution itself).

import numpy as np
import pandas as pd
from scipy import stats

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any()

# PINNED: the groups are split by arm LABEL, never passed as a numeric column to a regression, for
# the same reason r.R pins factor(arm).
groups = [d.y[d.arm == a].to_numpy() for a in (1, 2, 3)]
anova = stats.f_oneway(*groups)
tukey = stats.tukey_hsd(*groups)
ci = tukey.confidence_interval(confidence_level=0.95)

n, k = len(d), len(groups)
pooled_var = sum(((g - g.mean()) ** 2).sum() for g in groups) / (n - k)

print(f"F = {anova.statistic:.3f} on {k - 1} and {n - k} df")
print(tukey)

print("\n--- HARNESS ---")
print(f"f_statistic={anova.statistic:.10f}")
print(f"df_between={k - 1}\ndf_within={n - k}")
print(f"pooled_sd={float(np.sqrt(pooled_var)):.10f}")
# statistic[i, j] is mean(group i) - mean(group j), so [1, 0] is arm 2 minus arm 1.
for (i, j) in ((1, 0), (2, 0), (2, 1)):
    key = f"diff_{i + 1}_{j + 1}"
    print(f"{key}={tukey.statistic[i, j]:.10f}")
    print(f"{key}_lcl={ci.low[i, j]:.10f}")
    print(f"{key}_ucl={ci.high[i, j]:.10f}")
print(f"n={n}")
