# Cox regression with a time-varying exposure, in counting-process (start, stop] form
#
# The Python equivalent of r.R, using lifelines' CoxTimeVaryingFitter (Efron ties), independent of
# R's coxph.

import numpy as np
import pandas as pd
from lifelines import CoxTimeVaryingFitter
from scipy import stats

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any() and (d["stop"] > d["start"]).all()

# PINNED: the start/stop form; a CoxPHFitter on one row per person could only use baseline exposure.
fit = CoxTimeVaryingFitter().fit(d, id_col="id", event_col="event", start_col="start", stop_col="stop")

z = stats.norm.ppf(0.975)
b, se = fit.params_["exposed"], fit.standard_errors_["exposed"]
print(f"hazard ratio for current exposure {np.exp(b):.3f}, 95% {np.exp(b - z * se):.3f} to {np.exp(b + z * se):.3f}")

print("\n--- HARNESS ---")
print(f"log_hazard_ratio={b:.10f}\nlog_hazard_ratio_se={se:.10f}")
print(f"log_hazard_ratio_lcl={b - z * se:.10f}\nlog_hazard_ratio_ucl={b + z * se:.10f}")
print(f"events={int(d.event.sum())}\npeople={d.id.nunique()}")
