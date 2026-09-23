"""A prevalence from a cluster sample, with a design-based 95% interval.

An explicit with-replacement Taylor-linearized variance for a one-stage, equal-probability cluster
sample: the Python counterpart of R's survey package, written out rather than imported because no
standard Python package implements design-based survey variance. It is NOT a general survey
package: no strata, no finite population correction, no multistage sampling, no missing items.
"""

import numpy as np
import pandas as pd
from scipy import stats

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any() and (d.weight > 0).all()
assert set(d.outcome).issubset({0, 1})

# The patient-level prevalence is a RATIO of weighted totals, not the average of the clinic
# percentages. The two differ whenever clinic size is related to prevalence.
total = np.sum(d.weight)
prevalence = float(np.sum(d.weight * d.outcome) / total)

# Linearize, sum within each cluster, and take the between-cluster variance of those totals.
# PINNED: m/(m-1), the with-replacement factor R's survey package applies. Dropping it, or using
# the patient-level variance p(1-p)/n, gives a smaller and wrong standard error.
linearized = d.weight * (d.outcome - prevalence) / total
cluster_totals = d.assign(u=linearized).groupby("cluster", sort=True).u.sum()
m = len(cluster_totals)
variance = m / (m - 1) * np.sum((cluster_totals - cluster_totals.mean()) ** 2)
se = float(np.sqrt(variance))

# PINNED: a t reference with m - 1 degrees of freedom, where m is the number of clusters.
t = stats.t.ppf(0.975, m - 1)
lcl, ucl = prevalence - t * se, prevalence + t * se

print(f"{len(d)} patients in {m} clusters")
print(f"patient-level prevalence {prevalence:.4f}, design-based SE {se:.4f}, "
      f"95% {lcl:.4f} to {ucl:.4f}")

print("\n--- HARNESS ---")
print(f"prevalence={prevalence:.10f}")
print(f"prevalence_se={se:.10f}")
print(f"prevalence_lcl={lcl:.10f}")
print(f"prevalence_ucl={ucl:.10f}")
print(f"clusters={m}")
print(f"n={len(d)}")
