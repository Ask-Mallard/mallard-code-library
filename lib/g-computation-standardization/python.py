"""G-computation (marginal standardization): a population risk difference and risk ratio.

The same analysis as r.R: a logistic outcome model, both risks predicted for every patient and averaged,
and a standard error from the stacked estimating equations (model score plus the two standardized risks),
which includes the sampling of the covariates. The delta method on the model covariance alone, which
treats the covariates as fixed, is reported beside it. See r.R for why the coefficient is not the answer
and why the odds ratio is non-collapsible.
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

d = pd.read_csv("fixture.csv")
fit = smf.glm("event ~ treated + severity + comorbid", d, family=sm.families.Binomial()).fit()

d1, d0 = d.assign(treated=1), d.assign(treated=0)
p1 = np.asarray(fit.predict(d1))
p0 = np.asarray(fit.predict(d0))
mu1, mu0 = p1.mean(), p0.mean()
rd = mu1 - mu0

# PINNED: the stacked sandwich over (outcome-model coefficients, mu1, mu0).
cols = ["treated", "severity", "comorbid"]
design = lambda frame: np.column_stack([np.ones(len(frame))] + [frame[c].to_numpy(float) for c in cols])
X, X1, X0 = design(d), design(d1), design(d0)
p = np.asarray(fit.fittedvalues)
yv = d["event"].to_numpy(float)
n, k = X.shape
psi = np.column_stack([X * (yv - p)[:, None], p1 - mu1, p0 - mu0])
A = np.zeros((k + 2, k + 2))
A[:k, :k] = -(X * (p * (1 - p))[:, None]).T @ X
A[k, :k] = (X1 * (p1 * (1 - p1))[:, None]).sum(axis=0)
A[k, k] = -n
A[k + 1, :k] = (X0 * (p0 * (1 - p0))[:, None]).sum(axis=0)
A[k + 1, k + 1] = -n
Ainv = np.linalg.inv(A)
V = Ainv @ (psi.T @ psi) @ Ainv.T
g_rd = np.r_[np.zeros(k), 1.0, -1.0]
g_lrr = np.r_[np.zeros(k), 1 / mu1, -1 / mu0]
rd_se = float(np.sqrt(g_rd @ V @ g_rd))
lrr = float(np.log(mu1 / mu0))
lrr_se = float(np.sqrt(g_lrr @ V @ g_lrr))

grad = (X1 * (p1 * (1 - p1))[:, None]).mean(axis=0) - (X0 * (p0 * (1 - p0))[:, None]).mean(axis=0)
info_inv = np.linalg.inv((X * (p * (1 - p))[:, None]).T @ X)  # the model covariance, written out
rd_se_conditional = float(np.sqrt(grad @ info_inv @ grad))

crude = d.loc[d.treated == 1, "event"].mean() - d.loc[d.treated == 0, "event"].mean()
z = stats.norm.ppf(0.975)
print(f"risk difference {rd:.4f} (95% CI {rd - z * rd_se:.4f} to {rd + z * rd_se:.4f}); risk ratio {np.exp(lrr):.3f}")

print("\n--- HARNESS ---")
print(f"risk_treated={mu1:.10f}\nrisk_untreated={mu0:.10f}")
print(f"risk_difference={rd:.10f}\nrisk_difference_se={rd_se:.10f}")
print(f"risk_difference_lcl={rd - z * rd_se:.10f}\nrisk_difference_ucl={rd + z * rd_se:.10f}")
print(f"risk_difference_se_conditional={rd_se_conditional:.10f}")
print(f"log_risk_ratio={lrr:.10f}\nlog_risk_ratio_se={lrr_se:.10f}")
print(f"log_marginal_odds_ratio={np.log(mu1 / (1 - mu1)) - np.log(mu0 / (1 - mu0)):.10f}")
print(f"log_conditional_odds_ratio={fit.params['treated']:.10f}")
print(f"crude_risk_difference={crude:.10f}\nn={n}")
