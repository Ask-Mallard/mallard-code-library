"""Interrupted time series: segmented regression with autocorrelation-robust (Newey-West) errors.

The same analysis as r.R: rate ~ month + post + months_since by OLS, with a HAC covariance at lag 3,
Bartlett weights and no small-sample correction (use_correction=False), which is what
sandwich::NeweyWest(lag = 3, prewhite = FALSE, adjust = FALSE) computes. See r.R for the two effects and
for why the lag is stated rather than chosen from the data.
"""

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

d = pd.read_csv("fixture.csv")
model = smf.ols("rate ~ month + post + months_since", d)
fit = model.fit(cov_type="HAC", cov_kwds={"maxlags": 3, "kernel": "bartlett", "use_correction": False})
ols = model.fit()
b, se = fit.params, fit.bse
r = ols.resid.to_numpy()
rho1 = float(np.sum(r[1:] * r[:-1]) / np.sum(r ** 2))

z = stats.norm.ppf(0.975)
print(f"level change {b['post']:.2f} (95% CI {b['post'] - z * se['post']:.2f} to {b['post'] + z * se['post']:.2f}); "
      f"slope change {b['months_since']:.3f} a month")

print("\n--- HARNESS ---")
print(f"level_change={b['post']:.10f}\nlevel_change_se={se['post']:.10f}")
print(f"slope_change={b['months_since']:.10f}\nslope_change_se={se['months_since']:.10f}")
print(f"pre_slope={b['month']:.10f}")
print(f"level_change_se_ols={ols.bse['post']:.10f}\nslope_change_se_ols={ols.bse['months_since']:.10f}")
print(f"residual_lag1_autocorrelation={rho1:.10f}\nn={len(d)}")
