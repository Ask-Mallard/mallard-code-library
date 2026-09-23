"""Developing a penalized prediction model and correcting its apparent performance by the bootstrap.

The same procedure as r.R, step for step: ridge logistic regression on standardized predictors
(intercept unpenalized, penalty lambda / 2 x sum(beta^2) on the log-likelihood), fitted by Newton-Raphson
to 1e-12; the penalty chosen by 5-fold cross-validated deviance over the same 12-value grid; Harrell's
bootstrap optimism correction with every step, including the penalty choice, repeated in each of the
fixture's 50 resamples; and the external validation cohort for comparison. scikit-learn's
LogisticRegression uses C = 1 / lambda on a sum-of-losses objective and its own solvers, so the fit is
written out. See r.R for the reasoning.
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm

d = pd.read_csv("fixture.csv")
dev, val = d[d.cohort == "development"].reset_index(drop=True), d[d.cohort == "validation"].reset_index(drop=True)
xs = [f"x{j}" for j in range(1, 11)]
boot = [np.array(list(map(int, s.split()))) for s in pd.read_csv("bootstrap.csv")["rows"]]
lambdas = np.exp(np.linspace(np.log(0.1), np.log(100), 12))


def ridge(X, y, lam):
    mu, s = X.mean(axis=0), X.std(axis=0, ddof=1)
    Z = np.column_stack([np.ones(len(X)), (X - mu) / s])
    pen = np.diag(np.r_[0.0, np.full(X.shape[1], lam)])
    b = np.zeros(Z.shape[1])
    while True:
        p = 1 / (1 + np.exp(-Z @ b))
        step = np.linalg.solve((Z * (p * (1 - p))[:, None]).T @ Z + pen, Z.T @ (y - p) - pen @ b)
        b = b + step
        if np.max(np.abs(step)) < 1e-12:
            break
    beta = b[1:] / s
    return np.r_[b[0] - np.sum(beta * mu), beta]


lp = lambda coef, X: coef[0] + X @ coef[1:]


def auc(eta, y):
    x, z = eta[y == 1], eta[y == 0]
    return ((x[:, None] > z[None, :]) + 0.5 * (x[:, None] == z[None, :])).mean()


def slope(eta, y):
    return sm.GLM(y, sm.add_constant(eta), family=sm.families.Binomial()).fit(tol=1e-14, maxiter=100).params[1]


def deviance(eta, y):
    return -2 * np.sum(y * eta - np.log1p(np.exp(eta)))


def choose_lambda(X, y, fold):
    cv = [sum(deviance(lp(ridge(X[fold != f], y[fold != f], lam), X[fold == f]), y[fold == f])
              for f in np.unique(fold)) for lam in lambdas]
    return lambdas[int(np.argmin(cv))]


X, y = dev[xs].to_numpy(float), dev["event"].to_numpy(float)
lam = choose_lambda(X, y, dev["fold"].to_numpy())
fit = ridge(X, y, lam)
app_auc, app_slope = auc(lp(fit, X), y), slope(lp(fit, X), y)

opt = []
for idx in boot:
    Xb, yb = X[idx], y[idx]
    fb = ridge(Xb, yb, choose_lambda(Xb, yb, np.arange(len(idx)) % 5 + 1))
    opt.append((auc(lp(fb, Xb), yb) - auc(lp(fb, X), y), slope(lp(fb, Xb), yb) - slope(lp(fb, X), y)))
opt = np.array(opt)
cor_auc, cor_slope = app_auc - opt[:, 0].mean(), app_slope - opt[:, 1].mean()
Xv, yv = val[xs].to_numpy(float), val["event"].to_numpy(float)
ext_auc, ext_slope = auc(lp(fit, Xv), yv), slope(lp(fit, Xv), yv)

print(f"lambda {lam:.3f}; AUC apparent {app_auc:.3f}, corrected {cor_auc:.3f}, external {ext_auc:.3f}; "
      f"slope apparent {app_slope:.3f}, corrected {cor_slope:.3f}, external {ext_slope:.3f}")

print("\n--- HARNESS ---")
print(f"lambda={lam:.10f}")
print(f"apparent_auc={app_auc:.10f}\noptimism_auc={opt[:, 0].mean():.10f}\ncorrected_auc={cor_auc:.10f}\nexternal_auc={ext_auc:.10f}")
print(f"apparent_slope={app_slope:.10f}\noptimism_slope={opt[:, 1].mean():.10f}\ncorrected_slope={cor_slope:.10f}\nexternal_slope={ext_slope:.10f}")
print(f"corrected_minus_external_auc={cor_auc - ext_auc:.10f}\ncorrected_minus_external_slope={cor_slope - ext_slope:.10f}")
print(f"coef_x1={fit[1]:.10f}\nn_development={len(dev)}")
