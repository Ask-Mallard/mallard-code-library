"""Restricted mean survival time (RMST) and the difference between arms, up to a stated horizon.

The RMST is the area under the Kaplan-Meier curve from 0 to tau, computed here from lifelines' KM fit.
Its standard error is NOT taken from lifelines: restricted_mean_survival_time(..., return_variance=True)
returns the variance of a patient's restricted survival time, min(T, tau), the spread of the data, not
the sampling variance of the RMST estimate. Used as a standard error it gave 1.24 years where the
estimate's standard error is 0.072. The estimate's variance is computed here with the Greenwood-type
formula (Klein and Moeschberger, section 4.5): sum over event times t_j <= tau of
A_j^2 x d_j / (n_j (n_j - d_j)), where A_j is the area under the curve from t_j to tau.
"""

import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.utils import restricted_mean_survival_time
from scipy import stats

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any()


def rmst_with_se(time, event, tau):
    kmf = KaplanMeierFitter().fit(time, event)
    area = restricted_mean_survival_time(kmf, t=tau)
    # Step function of the KM estimate at the distinct event times up to tau.
    table = kmf.event_table
    times = table.index.to_numpy(float)
    events = table["observed"].to_numpy(float)
    at_risk = table["at_risk"].to_numpy(float)
    keep = (events > 0) & (times <= tau)
    t_j, d_j, n_j = times[keep], events[keep], at_risk[keep]
    surv = kmf.survival_function_at_times(t_j).to_numpy(float)          # S(t_j), just after t_j
    # Area from each event time to tau: the curve is flat at S(t_j) until the next event time.
    edges = np.append(t_j, tau)
    pieces = surv * np.diff(edges)
    area_after = np.cumsum(pieces[::-1])[::-1]                            # A_j
    var = np.sum(area_after ** 2 * d_j / (n_j * (n_j - d_j)))
    return float(area), float(np.sqrt(var))


# PINNED: tau = 4 years, stated in advance, and the estimator's own variance.
tau = 4.0
r1, s1 = rmst_with_se(d.time[d.arm == 1], d.event[d.arm == 1], tau)
r0, s0 = rmst_with_se(d.time[d.arm == 0], d.event[d.arm == 0], tau)
diff, se = r1 - r0, float(np.sqrt(s1 ** 2 + s0 ** 2))
z = stats.norm.ppf(0.975)
print(f"RMST to {tau:g} years: treated {r1:.3f}, control {r0:.3f}; difference {diff:.3f}, "
      f"95% {diff - z * se:.3f} to {diff + z * se:.3f}")

print("\n--- HARNESS ---")
print(f"rmst_treated={r1:.10f}\nrmst_treated_se={s1:.10f}")
print(f"rmst_control={r0:.10f}\nrmst_control_se={s0:.10f}")
print(f"rmst_difference={diff:.10f}\nrmst_difference_se={se:.10f}")
print(f"rmst_difference_lcl={diff - z * se:.10f}\nrmst_difference_ucl={diff + z * se:.10f}")
print(f"n={len(d)}")
