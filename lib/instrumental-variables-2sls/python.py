"""Instrumental variables: two-stage least squares (Mendelian-randomization style).

The same analysis as r.R, written out: the first stage projects the exposure on the instrument, the
second regresses the outcome on the projection, and the RESIDUALS for the standard error are the outcome
minus the exposure (not its projection) times the coefficient, as ivreg computes them. Running the two
stages as separate OLS fits gets the coefficient right and the standard error wrong. Pinned: the HC1
sandwich; the conventional SE and the first-stage F are reported beside it. See r.R for the three
assumptions.
"""

import numpy as np
import pandas as pd
from scipy import stats

d = pd.read_csv("fixture.csv")
n = len(d)
X = np.column_stack([np.ones(n), d["ldl"]])
Z = np.column_stack([np.ones(n), d["allele_score"]])
y = d["y"].to_numpy(float)

# First stage and projection.
Xhat = Z @ np.linalg.solve(Z.T @ Z, Z.T @ X)
beta = np.linalg.solve(Xhat.T @ X, Xhat.T @ y)
resid = y - X @ beta  # PINNED: structural residuals, with the observed exposure
k = X.shape[1]

bread = np.linalg.inv(Xhat.T @ Xhat)
meat = (Xhat * resid[:, None] ** 2).T @ Xhat
V_hc1 = bread @ meat @ bread * n / (n - k)
V_conv = bread * np.sum(resid ** 2) / (n - k)
se_hc1, se_conv = float(np.sqrt(V_hc1[1, 1])), float(np.sqrt(V_conv[1, 1]))

# First-stage F for the single instrument (conventional), as ivreg's weak-instrument diagnostic.
g = np.linalg.solve(Z.T @ Z, Z.T @ X[:, 1])
r1 = X[:, 1] - Z @ g
V1 = np.linalg.inv(Z.T @ Z) * np.sum(r1 ** 2) / (n - 2)
first_f = float(g[1] ** 2 / V1[1, 1])

ols = float(np.linalg.lstsq(X, y, rcond=None)[0][1])
b = float(beta[1])
z = stats.norm.ppf(0.975)
print(f"2SLS effect {b:.3f} (95% CI {b - z * se_hc1:.3f} to {b + z * se_hc1:.3f}); OLS {ols:.3f}; first-stage F {first_f:.1f}")

print("\n--- HARNESS ---")
print(f"iv_effect={b:.10f}\niv_se={se_hc1:.10f}\niv_lcl={b - z * se_hc1:.10f}\niv_ucl={b + z * se_hc1:.10f}")
print(f"iv_se_conventional={se_conv:.10f}")
print(f"first_stage_f={first_f:.10f}\nols_effect={ols:.10f}\nn={n}")
