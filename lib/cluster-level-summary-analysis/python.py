# Cluster-level summary analysis for a cluster trial with few clusters
#
# The Python equivalent of r.R, using scipy's two-sample t-test on the cluster proportions.

import pandas as pd
from scipy import stats

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any()

# PINNED: one value per CLUSTER; equal_var=True (the classical t on k - 2 df), intervention first.
clusters = d.groupby(["cluster", "arm"], as_index=False).infection.mean()
treated = clusters.infection[clusters.arm == 1]
control = clusters.infection[clusters.arm == 0]
fit = stats.ttest_ind(treated, control, equal_var=True)
ci = fit.confidence_interval(confidence_level=0.95)
diff = treated.mean() - control.mean()

print(f"{len(clusters)} clusters; difference in mean cluster proportion {diff:.4f}, "
      f"95% {ci.low:.4f} to {ci.high:.4f}, {int(fit.df)} df")

print("\n--- HARNESS ---")
print(f"difference={diff:.10f}\ndifference_lcl={ci.low:.10f}\ndifference_ucl={ci.high:.10f}")
print(f"df={int(fit.df)}\nclusters={len(clusters)}\nn={len(d)}")
