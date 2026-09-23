"""Controlled interrupted time series: an intervention site against a control site.

The same analysis as r.R: segmented regression of the difference series (intervention minus control)
with Newey-West errors at lag 3, Bartlett weights, no small-sample correction. The single-series analysis
of the intervention site is reported beside it. See r.R for why the difference series gives the full
interaction model's coefficients and what the control series assumes.
"""

import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

d = pd.read_csv("fixture.csv")
d["difference"] = d["intervention_rate"] - d["control_rate"]
fit = smf.ols("difference ~ month + post + months_since", d).fit(
    cov_type="HAC", cov_kwds={"maxlags": 3, "kernel": "bartlett", "use_correction": False})
b, se = fit.params, fit.bse
single = smf.ols("intervention_rate ~ month + post + months_since", d).fit().params

z = stats.norm.ppf(0.975)
print(f"level change relative to control {b['post']:.2f} (95% CI {b['post'] - z * se['post']:.2f} to "
      f"{b['post'] + z * se['post']:.2f}); slope change {b['months_since']:.3f}")

print("\n--- HARNESS ---")
print(f"level_change={b['post']:.10f}\nlevel_change_se={se['post']:.10f}")
print(f"slope_change={b['months_since']:.10f}\nslope_change_se={se['months_since']:.10f}")
print(f"single_series_level_change={single['post']:.10f}\nsingle_series_slope_change={single['months_since']:.10f}")
print(f"n={len(d)}")
