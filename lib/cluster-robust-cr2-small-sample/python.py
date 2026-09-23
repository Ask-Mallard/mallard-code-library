"""Small-sample cluster-robust inference: CR2 standard errors with Satterthwaite degrees of freedom.

No standard Python package implements CR2 with Satterthwaite degrees of freedom (statsmodels' cluster
option is the CR1-style sandwich with G - 1 df). Both are written out here for an ordinary least-squares
fit, following Bell and McCaffrey (2002) and Pustejovsky and Tipton (2018):

  CR2:  V = (X'X)^-1 [ sum_g X_g' A_g e_g e_g' A_g X_g ] (X'X)^-1,  A_g = (I - H_gg)^(-1/2)
  df:   for the contrast c, with a_g = A_g X_g (X'X)^-1 c placed in cluster g's rows and B = sum_g a_g a_g',
        df = tr(B M)^2 / tr(B M B M),  M = I - H  (the Satterthwaite approximation under a working
        model of independent, equal-variance errors)
"""

import numpy as np
import pandas as pd
from scipy import stats

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any()

X = np.column_stack([np.ones(len(d)), d.treated.to_numpy(float), (d.age - 50).to_numpy(float)])
y = d.y.to_numpy(float)
XtX_inv = np.linalg.inv(X.T @ X)
beta = XtX_inv @ X.T @ y
e = y - X @ beta
H = X @ XtX_inv @ X.T
M = np.eye(len(d)) - H


def inv_sqrt(S):
    vals, vecs = np.linalg.eigh(S)
    return vecs @ np.diag(1 / np.sqrt(vals)) @ vecs.T


# PINNED: CR2, the bias-reduced adjustment, not the plain sandwich.
c = np.array([0.0, 1.0, 0.0])                   # the treatment coefficient
meat = np.zeros((3, 3))
B = np.zeros((len(d), len(d)))
for idx in d.groupby("clinic").indices.values():
    A = inv_sqrt(np.eye(len(idx)) - H[np.ix_(idx, idx)])
    u = X[idx].T @ A @ e[idx]
    meat += np.outer(u, u)
    a = np.zeros(len(d))
    a[idx] = A @ X[idx] @ XtX_inv @ c
    B += np.outer(a, a)
V = XtX_inv @ meat @ XtX_inv
se = float(np.sqrt(V[1, 1]))
BM = B @ M
df = float(np.trace(BM) ** 2 / np.trace(BM @ BM))

tq = stats.t.ppf(0.975, df)
print(f"treatment effect {beta[1]:.3f}, CR2 SE {se:.4f}, Satterthwaite df {df:.2f}, "
      f"95% {beta[1] - tq * se:.3f} to {beta[1] + tq * se:.3f}")

print("\n--- HARNESS ---")
print(f"effect={beta[1]:.10f}\neffect_cr2_se={se:.10f}\neffect_df={df:.10f}")
print(f"effect_lcl={beta[1] - tq * se:.10f}\neffect_ucl={beta[1] + tq * se:.10f}")
print(f"clusters={d.clinic.nunique()}\nn={len(d)}")
