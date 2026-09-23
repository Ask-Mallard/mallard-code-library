# Weibull regression: an accelerated failure time (time ratio), and the equivalent hazard ratio
#
# The Python equivalent of r.R, using lifelines' WeibullAFTFitter, independent of R's survreg.

import numpy as np
import pandas as pd
from lifelines import WeibullAFTFitter
from scipy import stats

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any()

# lifelines stops about 5e-6 short of the maximum on this fixture (an independent optimisation agrees
# with survreg to 3e-8); its fit_options did not tighten that. The gap is inside the library's standard
# 1e-4 agreement tolerance, and far below anything clinical.
fit = WeibullAFTFitter().fit(d[["time", "event", "treated"]], duration_col="time", event_col="event")

# PINNED: lifelines' lambda_ coefficients are survreg's log time ratios; its rho_ is the Weibull SHAPE,
# modelled on the log scale, so shape = exp(rho_ intercept) and survreg's scale = 1 / shape.
b = fit.params_[("lambda_", "treated")]
se = fit.standard_errors_[("lambda_", "treated")]
shape = float(np.exp(fit.params_[("rho_", "Intercept")]))
z = stats.norm.ppf(0.975)
print(f"time ratio {np.exp(b):.3f}, 95% {np.exp(b - z * se):.3f} to {np.exp(b + z * se):.3f}; "
      f"shape {shape:.3f}; hazard ratio {np.exp(-b * shape):.3f}")

print("\n--- HARNESS ---")
print(f"log_time_ratio={b:.10f}\nlog_time_ratio_se={se:.10f}")
print(f"log_time_ratio_lcl={b - z * se:.10f}\nlog_time_ratio_ucl={b + z * se:.10f}")
print(f"shape={shape:.10f}")
print(f"log_hazard_ratio={-b * shape:.10f}")
print(f"n={len(d)}")
