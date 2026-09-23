# Paired t-test: the mean within-person change, after minus before, with a 95% interval
#
# The Python equivalent of r.R, using scipy's ttest_rel, an implementation independent of R's t.test.

import numpy as np
import pandas as pd
from scipy import stats

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any()

# PINNED: ttest_rel (paired), with after FIRST so the difference is after minus before.
# ttest_ind on the same columns would be the unpaired test, and a reversed order flips the sign.
fit = stats.ttest_rel(d.after, d.before)
ci = fit.confidence_interval(confidence_level=0.95)
change = d.after - d.before

print(f"{len(d)} pairs; mean change (after - before) {change.mean():.2f}, "
      f"95% {ci.low:.2f} to {ci.high:.2f}; t = {fit.statistic:.3f} on {int(fit.df)} df")

print("\n--- HARNESS ---")
print(f"mean_change={change.mean():.10f}")
print(f"mean_change_lcl={ci.low:.10f}")
print(f"mean_change_ucl={ci.high:.10f}")
print(f"sd_change={float(np.std(change, ddof=1)):.10f}")
print(f"t_statistic={fit.statistic:.10f}")
print(f"df={int(fit.df)}")
print(f"n={len(d)}")
