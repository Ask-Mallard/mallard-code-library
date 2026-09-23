"""Causal mediation with a binary outcome: natural direct and indirect effects by the mediation formula.

The same analysis as r.R: a logistic mediator model and a logistic outcome model WITH the
treatment-by-mediator interaction, the potential risks E[Y(a, M(a*))] computed by summing exactly over
the binary mediator and averaging over the covariate, and a standard error from the stacked estimating
equations. The difference-method log odds ratio is reported beside it, not as an indirect effect: see
r.R for why Baron-Kenny and product-of-coefficients are invalid with a binary outcome.
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
from scipy.special import expit, logit

d = pd.read_csv("fixture.csv")
n = len(d)
fit = lambda f: smf.glm(f, d, family=sm.families.Binomial()).fit(tol=1e-12)
mfit = fit("mediator ~ treated + c")
yfit = fit("event ~ treated * mediator + c")
g = mfit.params.to_numpy()
b = yfit.params.to_numpy()
assert list(yfit.params.index) == ["Intercept", "treated", "mediator", "treated:mediator", "c"]

c = d["c"].to_numpy(float)
one = np.ones(n)
XM = lambda a: np.column_stack([one, a * one, c])
XY = lambda a, m: np.column_stack([one, a * one, m * one, a * m * one, c])  # column order of yfit.params
pm = lambda a: expit(XM(a) @ g)
py = lambda a, m: expit(XY(a, m) @ b)
mu_i = lambda a, a_star: py(a, 1) * pm(a_star) + py(a, 0) * (1 - pm(a_star))

# PINNED: the mediation formula, summed exactly over the binary mediator.
pairs = [(1, 0), (0, 0), (1, 1)]
mu = np.array([mu_i(a, s).mean() for a, s in pairs])
nde, nie = mu[0] - mu[1], mu[2] - mu[0]

# PINNED: the stacked sandwich over (mediator coefficients, outcome coefficients, three risks).
Xm = np.column_stack([one, d["treated"], c])
Xy = XY(d["treated"].to_numpy(float), d["mediator"].to_numpy(float))
p_m, p_y = expit(Xm @ g), expit(Xy @ b)
k1, k2 = Xm.shape[1], Xy.shape[1]
psi = np.column_stack([Xm * (d["mediator"].to_numpy() - p_m)[:, None], Xy * (d["event"].to_numpy() - p_y)[:, None]]
                      + [mu_i(a, s) - mu[j] for j, (a, s) in enumerate(pairs)])
K = k1 + k2 + 3
A = np.zeros((K, K))
A[:k1, :k1] = -(Xm * (p_m * (1 - p_m))[:, None]).T @ Xm
A[k1:k1 + k2, k1:k1 + k2] = -(Xy * (p_y * (1 - p_y))[:, None]).T @ Xy
for j, (a, s) in enumerate(pairs):
    q, y1, y0 = pm(s), py(a, 1), py(a, 0)
    A[k1 + k2 + j, :k1] = (XM(s) * ((y1 - y0) * q * (1 - q))[:, None]).sum(axis=0)
    A[k1 + k2 + j, k1:k1 + k2] = (XY(a, 1) * (y1 * (1 - y1) * q)[:, None]
                                  + XY(a, 0) * (y0 * (1 - y0) * (1 - q))[:, None]).sum(axis=0)
    A[k1 + k2 + j, k1 + k2 + j] = -n
Ainv = np.linalg.inv(A)
V = Ainv @ (psi.T @ psi) @ Ainv.T
se = lambda w: float(np.sqrt(np.r_[np.zeros(k1 + k2), w] @ V @ np.r_[np.zeros(k1 + k2), w]))
nde_se, nie_se = se([1, -1, 0]), se([-1, 0, 1])

total_coef = fit("event ~ treated + c").params["treated"]
direct_coef = fit("event ~ treated + mediator + c").params["treated"]
nie_log_or = logit(mu[2]) - logit(mu[0])

z = stats.norm.ppf(0.975)
print(f"NDE {nde:.4f} (95% CI {nde - z * nde_se:.4f} to {nde + z * nde_se:.4f}); "
      f"NIE {nie:.4f} ({nie - z * nie_se:.4f} to {nie + z * nie_se:.4f}); proportion mediated {nie / (nde + nie):.2f}")

print("\n--- HARNESS ---")
print(f"natural_direct_effect={nde:.10f}\nnatural_direct_effect_se={nde_se:.10f}")
print(f"natural_indirect_effect={nie:.10f}\nnatural_indirect_effect_se={nie_se:.10f}")
print(f"total_effect={nde + nie:.10f}\ntotal_effect_se={se([0, -1, 1]):.10f}")
print(f"nie_log_odds_ratio={nie_log_or:.10f}\ndifference_method_log_odds_ratio={total_coef - direct_coef:.10f}")
print(f"n={n}")
