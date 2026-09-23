"""Median (quantile) regression: the adjusted difference in median length of stay.

statsmodels' QuantReg fits by iteratively reweighted least squares, an approximation to the linear
program, and its default standard error is a kernel sandwich unlike any of R's. So the fit here is
the exact linear program (scipy's HiGHS), and the standard error is the Hendricks-Koenker local
sandwich with the Hall-Sheather bandwidth, written out following Koenker (2005) section 3.4, the
estimator R's summary.rq(se = "nid") reports.
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import linprog

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any()

X = np.column_stack([np.ones(len(d)), d.treated.to_numpy(float), (d.age - 60).to_numpy(float)])
y = d.los_days.to_numpy(float)
n, p = X.shape


def rq_fit(tau):
    """Exact quantile regression: minimise sum tau*u+ + (1-tau)*u- subject to X b + u+ - u- = y."""
    c = np.concatenate([np.zeros(p), np.full(n, tau), np.full(n, 1 - tau)])
    A = np.hstack([X, np.eye(n), -np.eye(n)])
    bounds = [(None, None)] * p + [(0, None)] * (2 * n)
    res = linprog(c, A_eq=A, b_eq=y, bounds=bounds, method="highs")
    assert res.status == 0, res.message
    return res.x[:p]


tau = 0.5
beta = rq_fit(tau)

# Hall-Sheather bandwidth (alpha = 0.05), then the sparsity from fits at tau +/- h.
x0 = stats.norm.ppf(tau)
h = n ** (-1 / 3) * stats.norm.ppf(0.975) ** (2 / 3) * ((1.5 * stats.norm.pdf(x0) ** 2) / (2 * x0 ** 2 + 1)) ** (1 / 3)
dyhat = X @ (rq_fit(tau + h) - rq_fit(tau - h))
f = np.maximum(0, 2 * h / (dyhat - np.sqrt(np.finfo(float).eps)))
fxx_inv = np.linalg.inv((X * f[:, None]).T @ X)
cov = tau * (1 - tau) * fxx_inv @ (X.T @ X) @ fxx_inv
se = np.sqrt(np.diag(cov))

tq = stats.t.ppf(0.975, n - p)
print(np.column_stack([beta, se]))

print("\n--- HARNESS ---")
print(f"median_effect={beta[1]:.10f}\nmedian_effect_se={se[1]:.10f}")
print(f"median_effect_lcl={beta[1] - tq * se[1]:.10f}\nmedian_effect_ucl={beta[1] + tq * se[1]:.10f}")
print(f"age_slope={beta[2]:.10f}")
print(f"n={n}")
