"""Overlap weights: the average treatment effect in the overlap population (ATO).

The same analysis as r.R, point for point: weights 1 - e for the treated and e for the controls, the
weighted difference in means, and a standard error from the stacked estimating equations (the propensity
model's score and the two weighted means), so the estimation of the score is accounted for. The HC0
sandwich that treats the weights as known is reported beside it. See r.R for why the target is the ATO
and why the weighted covariate means balance exactly.
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

d = pd.read_csv("fixture.csv")
a = d["treated"].to_numpy(float)
y = d["y"].to_numpy(float)
X = np.column_stack([np.ones(len(d)), d["severity"], d["comorbid"]])

ps = smf.glm("treated ~ severity + comorbid", d, family=sm.families.Binomial()).fit()
e = np.asarray(ps.fittedvalues)
w = np.where(a == 1, 1 - e, e)  # PINNED: overlap weights

mu1 = np.sum(w * a * y) / np.sum(w * a)
mu0 = np.sum(w * (1 - a) * y) / np.sum(w * (1 - a))
ato = mu1 - mu0

# PINNED: the stacked M-estimation sandwich over (propensity coefficients, mu1, mu0).
n, k = X.shape
psi = np.column_stack([X * (a - e)[:, None], w * a * (y - mu1), w * (1 - a) * (y - mu0)])
de = e * (1 - e)
A = np.zeros((k + 2, k + 2))
A[:k, :k] = -(X * de[:, None]).T @ X
A[k, :k] = np.sum((-a * (y - mu1) * de)[:, None] * X, axis=0)
A[k, k] = -np.sum(w * a)
A[k + 1, :k] = np.sum(((1 - a) * (y - mu0) * de)[:, None] * X, axis=0)
A[k + 1, k + 1] = -np.sum(w * (1 - a))
Ainv = np.linalg.inv(A)
V = Ainv @ (psi.T @ psi) @ Ainv.T
contrast = np.r_[np.zeros(k), 1.0, -1.0]
se = float(np.sqrt(contrast @ V @ contrast))

se_fixed = float(np.sqrt(np.sum((w * a * (y - mu1)) ** 2) / np.sum(w * a) ** 2
                         + np.sum((w * (1 - a) * (y - mu0)) ** 2) / np.sum(w * (1 - a)) ** 2))

sev = d["severity"].to_numpy()
wmean = lambda x, g: np.sum(w * g * x) / np.sum(w * g)
smd_w = (wmean(sev, a) - wmean(sev, 1 - a)) / np.sqrt((sev[a == 1].var(ddof=1) + sev[a == 0].var(ddof=1)) / 2)

z = stats.norm.ppf(0.975)
print(f"ATO {ato:.3f} (95% CI {ato - z * se:.3f} to {ato + z * se:.3f})")

print("\n--- HARNESS ---")
print(f"ato={ato:.10f}\nato_se={se:.10f}\nato_lcl={ato - z * se:.10f}\nato_ucl={ato + z * se:.10f}")
print(f"ato_se_weights_fixed={se_fixed:.10f}")
print(f"smd_severity_weighted={smd_w:.10f}")
print(f"max_weight={w.max():.10f}\nn={n}")
