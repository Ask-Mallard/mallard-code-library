# Cause-specific Cox regression under competing risks
#
# The Python equivalent of r.R, using lifelines' CoxPHFitter (Efron ties), independent of R's coxph.

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from scipy import stats

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any() and set(d.status).issubset({0, 1, 2})

# PINNED: the event indicator is status == 1; deaths (status 2) are censored for THIS model only.
cs = pd.DataFrame({"time": d.time, "relapse": (d.status == 1).astype(int), "treated": d.treated})
fit = CoxPHFitter().fit(cs, duration_col="time", event_col="relapse")

z = stats.norm.ppf(0.975)
b, se = fit.params_["treated"], fit.standard_errors_["treated"]
print(f"cause-specific HR for relapse {np.exp(b):.3f}, 95% {np.exp(b - z * se):.3f} to {np.exp(b + z * se):.3f}")

print("\n--- HARNESS ---")
print(f"log_cs_hazard_ratio={b:.10f}\nlog_cs_hazard_ratio_se={se:.10f}")
print(f"log_cs_hazard_ratio_lcl={b - z * se:.10f}\nlog_cs_hazard_ratio_ucl={b + z * se:.10f}")
print(f"relapses={int((d.status == 1).sum())}\ncompeting_deaths={int((d.status == 2).sum())}\nn={len(d)}")
