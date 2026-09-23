"""Self-controlled case series (SCCS): conditional Poisson regression within each case.

The same analysis as r.R: a Poisson regression with one fixed effect per case and log(interval length)
as the offset, which has the same estimate and standard error as the SCCS conditional likelihood. Age
band is in the model; the fit without it is reported beside it. See r.R for the design's assumptions.

The fit is iteratively reweighted least squares WRITTEN OUT, each step a linear solve of the normal
equations, to a deviance change below 1e-12. statsmodels' GLM was used first: its IRLS solves each step
by SVD, which failed to converge ("SVD did not converge") on 1 of the 100 calibration fixtures locally,
and with wls_method="qr" its final step still calls pinv (an SVD), which failed intermittently on CI's
Linux runner for fixtures it had fitted before. A linear solve has no such failure mode here.
"""

import numpy as np
import pandas as pd
from scipy import stats

d = pd.read_csv("fixture.csv")
off = np.log(d["days"].to_numpy(float))
y = d["events"].to_numpy(float)
children = pd.get_dummies(d["child"], prefix="child", drop_first=True, dtype=float)
ages = pd.get_dummies(d["age_band"], prefix="age", drop_first=True, dtype=float)


def poisson_irls(X):
    """Poisson log-link IRLS with an offset; returns coefficients and their model-based covariance."""
    beta = np.zeros(X.shape[1])
    beta[0] = np.log(y.sum() / np.exp(off).sum())
    dev_old = np.inf
    for _ in range(100):
        eta = X @ beta + off
        mu = np.exp(eta)
        z = eta - off + (y - mu) / mu
        XtW = X.T * mu
        beta = np.linalg.solve(XtW @ X, XtW @ z)
        mu = np.exp(X @ beta + off)
        dev = 2 * np.sum(np.where(y > 0, y * np.log(np.where(y > 0, y, 1) / mu), 0) - (y - mu))
        if abs(dev - dev_old) < 1e-12 * (abs(dev) + 0.1):
            break
        dev_old = dev
    return beta, np.linalg.inv((X.T * mu) @ X)


X_full = np.column_stack([np.ones(len(d)), d["risk"], ages, children])
beta, cov = poisson_irls(X_full)
b, se = float(beta[1]), float(np.sqrt(cov[1, 1]))
no_age = float(poisson_irls(np.column_stack([np.ones(len(d)), d["risk"], children]))[0][1])

z = stats.norm.ppf(0.975)
print(f"incidence rate ratio {np.exp(b):.2f} (95% CI {np.exp(b - z * se):.2f} to {np.exp(b + z * se):.2f})")

print("\n--- HARNESS ---")
print(f"log_irr={b:.10f}\nlog_irr_se={se:.10f}\nlog_irr_lcl={b - z * se:.10f}\nlog_irr_ucl={b + z * se:.10f}")
print(f"log_irr_age_ignored={no_age:.10f}")
print(f"cases={d.child.nunique()}\nevents={int(d.events.sum())}")
