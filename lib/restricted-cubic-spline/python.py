"""Restricted (natural) cubic spline for a nonlinear continuous exposure.

The Python equivalent of r.R with a DIFFERENT basis: Harrell's truncated-power restricted cubic
spline, written out, where R uses the B-spline basis of splines::ns. The two bases span the same
functions for the same knots (cubic between the outer knots, linear beyond), so the coefficients differ
and every prediction and contrast must agree. That is why only predictions and contrasts are reported.
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any()

# PINNED: the same five knots as r.R, written out rather than taken from data quantiles.
KNOTS = np.array([30.0, 45.0, 60.0, 70.0, 80.0])


def rcs(x, t=KNOTS):
    x = np.asarray(x, dtype=float)
    k = len(t)
    pos = lambda v: np.maximum(v, 0.0) ** 3
    cols = [x]
    for j in range(k - 2):
        cols.append(pos(x - t[j]) - pos(x - t[k - 2]) * (t[k - 1] - t[j]) / (t[k - 1] - t[k - 2])
                    + pos(x - t[k - 1]) * (t[k - 2] - t[j]) / (t[k - 1] - t[k - 2]))
    return sm.add_constant(np.column_stack(cols), has_constant="add")


fit = sm.OLS(d.sbp.to_numpy(float), rcs(d.age)).fit()
Xn = rcs([40.0, 50.0, 70.0])
pred = Xn @ fit.params
contrast = Xn[2] - Xn[1]                     # age 70 minus age 50
est = float(contrast @ fit.params)
se = float(np.sqrt(contrast @ fit.cov_params() @ contrast))
tq = stats.t.ppf(0.975, fit.df_resid)

print(f"predicted SBP at 40, 50, 70: {pred[0]:.2f}, {pred[1]:.2f}, {pred[2]:.2f}")
print(f"age 70 vs 50: {est:.2f} mmHg, 95% {est - tq * se:.2f} to {est + tq * se:.2f}")

print("\n--- HARNESS ---")
print(f"pred_40={pred[0]:.10f}\npred_50={pred[1]:.10f}\npred_70={pred[2]:.10f}")
print(f"contrast_70_50={est:.10f}\ncontrast_70_50_se={se:.10f}")
print(f"contrast_70_50_lcl={est - tq * se:.10f}\ncontrast_70_50_ucl={est + tq * se:.10f}")
print(f"n={len(d)}")
