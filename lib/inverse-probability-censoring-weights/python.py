"""Inverse probability of censoring weights (IPCW) for informative dropout.

The same analysis as r.R: a logistic model for staying in the study given arm, the prognostic factor and
their interaction; completers weighted by 1 / P(stayed); each arm's weighted risk; and a standard error
from the stacked estimating equations (the dropout model's score and the two weighted risks). The SE that
treats the weights as known and the complete-case difference are reported beside it. See r.R for the
assumption this rests on.
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

d = pd.read_csv("fixture.csv")
d["stayed"] = d["event"].notna().astype(int)
cm = smf.glm("stayed ~ treated * l", d, family=sm.families.Binomial()).fit(tol=1e-12)
p = np.asarray(cm.fittedvalues)
X = np.column_stack([np.ones(len(d)), d["treated"], d["l"], d["treated"] * d["l"]])
a = d["treated"].to_numpy(float)
r = d["stayed"].to_numpy(float)
y = d["event"].fillna(0).to_numpy(float)
w = r / p  # PINNED: 1 / P(stayed | A, L) for completers, 0 for dropouts

mu1 = np.sum(w * a * y) / np.sum(w * a)
mu0 = np.sum(w * (1 - a) * y) / np.sum(w * (1 - a))
rd = mu1 - mu0

k = X.shape[1]
psi = np.column_stack([X * (r - p)[:, None], w * a * (y - mu1), w * (1 - a) * (y - mu0)])
A = np.zeros((k + 2, k + 2))
A[:k, :k] = -(X * (p * (1 - p))[:, None]).T @ X
dw = -(r * (1 - p) / p)
A[k, :k] = ((dw * a * (y - mu1))[:, None] * X).sum(axis=0)
A[k, k] = -np.sum(w * a)
A[k + 1, :k] = ((dw * (1 - a) * (y - mu0))[:, None] * X).sum(axis=0)
A[k + 1, k + 1] = -np.sum(w * (1 - a))
Ai = np.linalg.inv(A)
V = Ai @ (psi.T @ psi) @ Ai.T
c = np.r_[np.zeros(k), 1.0, -1.0]
se = float(np.sqrt(c @ V @ c))
se_fixed = float(np.sqrt(np.sum((w * a * (y - mu1)) ** 2) / np.sum(w * a) ** 2
                         + np.sum((w * (1 - a) * (y - mu0)) ** 2) / np.sum(w * (1 - a)) ** 2))

obs = r == 1
cc = y[obs & (a == 1)].mean() - y[obs & (a == 0)].mean()

z = stats.norm.ppf(0.975)
print(f"IPCW risk difference {rd:.4f} (95% CI {rd - z * se:.4f} to {rd + z * se:.4f}); complete-case {cc:.4f}")

print("\n--- HARNESS ---")
print(f"risk_difference={rd:.10f}\nrisk_difference_se={se:.10f}")
print(f"risk_difference_lcl={rd - z * se:.10f}\nrisk_difference_ucl={rd + z * se:.10f}")
print(f"risk_difference_se_weights_fixed={se_fixed:.10f}")
print(f"complete_case_risk_difference={cc:.10f}\nmax_weight={w.max():.10f}")
print(f"dropped={int((r == 0).sum())}\nn={len(d)}")
