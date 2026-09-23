"""Andersen-Gill model for recurrent events, with a robust (clustered) variance.

The rate ratio comes from lifelines' CoxPHFitter with delayed entry (entry_col). The robust variance
does NOT: in lifelines 0.28.0, CoxPHFitter(..., entry_col=..., cluster_col=..., robust=True) reported a
standard error of 0.0625 on this fixture, SMALLER than the naive 0.0725, where the Lin-Wei sandwich
clustered on the person is 0.0905 (R's coxph agrees). So the sandwich is written out here: each row's
score residual (its own event term minus its share of every risk set it sat in), summed by person, with
the information matrix from the partial likelihood. Exact for Efron when no event times tie, which is
asserted.
"""

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from scipy import stats

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any() and (d["stop"] > d["start"]).all()
assert not d.stop[d.event == 1].duplicated().any(), "tied event times: the written-out sandwich assumes none"

# PINNED: entry_col="start", the counting-process form. The point estimate is lifelines'.
fit = CoxPHFitter().fit(d[["treated", "start", "stop", "event"]], duration_col="stop", event_col="event",
                        entry_col="start")
b = float(fit.params_["treated"])
naive_se = float(fit.standard_errors_["treated"])

# PINNED: the robust variance, clustered on the person, written out (see the module docstring).
s, t, e = d.start.to_numpy(), d.stop.to_numpy(), d.event.to_numpy()
x = d.treated.to_numpy(float)
event_times = t[e == 1]
at_risk = (s[None, :] < event_times[:, None]) & (t[None, :] >= event_times[:, None])
risk = np.exp(b * x)
s0 = (at_risk * risk).sum(axis=1)
xbar = (at_risk * risk * x).sum(axis=1) / s0
information = ((at_risk * risk * x ** 2).sum(axis=1) / s0 - xbar ** 2).sum()
own = np.zeros(len(d))
own[e == 1] = x[e == 1] - xbar
share = (at_risk * risk * (x[None, :] - xbar[:, None]) / s0[:, None]).sum(axis=0)
per_person = pd.Series(own - share).groupby(d.id.to_numpy()).sum().to_numpy()
se = float(np.sqrt((per_person ** 2).sum()) / information)

z = stats.norm.ppf(0.975)
print(f"rate ratio {np.exp(b):.3f}, robust 95% {np.exp(b - z * se):.3f} to {np.exp(b + z * se):.3f} "
      f"(robust SE {se:.4f}, naive {naive_se:.4f})")

print("\n--- HARNESS ---")
print(f"log_rate_ratio={b:.10f}\nlog_rate_ratio_robust_se={se:.10f}")
print(f"log_rate_ratio_lcl={b - z * se:.10f}\nlog_rate_ratio_ucl={b + z * se:.10f}")
print(f"log_rate_ratio_naive_se={naive_se:.10f}")
print(f"events={int(d.event.sum())}\npeople={d.id.nunique()}")
