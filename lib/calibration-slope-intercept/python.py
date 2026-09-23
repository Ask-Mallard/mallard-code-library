"""Calibration of a published prediction model in new patients: slope, intercept, O/E.

The same analysis as r.R: the calibration slope from a logistic regression of the outcome on the
published linear predictor; calibration-in-the-large from a logistic regression with the linear
predictor as an OFFSET (slope fixed at 1), which is not the intercept of the slope model; the
observed-to-expected ratio; and the AUC for contrast. The calibration curve should be plotted (see
r.R) but is not compared across languages. The GLMs run to a tight tolerance.
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.special import expit

d = pd.read_csv("fixture.csv")
y = d["event"].to_numpy(float)
lp = d["lp"].to_numpy(float)
slope_fit = sm.GLM(y, sm.add_constant(lp), family=sm.families.Binomial()).fit(tol=1e-12)
citl_fit = sm.GLM(y, np.ones((len(y), 1)), family=sm.families.Binomial(), offset=lp).fit(tol=1e-12)
p = expit(lp)
oe = y.sum() / p.sum()
x, z = lp[y == 1], lp[y == 0]
auc = ((x[:, None] > z[None, :]) + 0.5 * (x[:, None] == z[None, :])).mean()

# The SEs written out as inverse information: statsmodels' IRLS covariance differed from R's vcov by 1e-6.
X = np.column_stack([np.ones(len(lp)), lp])
ps = expit(X @ slope_fit.params)
slope_se = float(np.sqrt(np.linalg.inv((X * (ps * (1 - ps))[:, None]).T @ X)[1, 1]))
pc = expit(citl_fit.params[0] + lp)
citl_se = float(1 / np.sqrt(np.sum(pc * (1 - pc))))

print(f"calibration slope {slope_fit.params[1]:.3f}, calibration-in-the-large {citl_fit.params[0]:.3f}, "
      f"O/E {oe:.3f}, AUC {auc:.3f}")

print("\n--- HARNESS ---")
print(f"calibration_slope={slope_fit.params[1]:.10f}\ncalibration_slope_se={slope_se:.10f}")
print(f"calibration_in_the_large={citl_fit.params[0]:.10f}\ncalibration_in_the_large_se={citl_se:.10f}")
print(f"slope_model_intercept={slope_fit.params[0]:.10f}")
print(f"observed_expected={oe:.10f}\nauc={auc:.10f}")
print(f"n={len(d)}\nevents={int(y.sum())}")
