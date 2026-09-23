# Zero-inflated Poisson regression for counts with structural zeros
#
# The Python equivalent of r.R, using statsmodels' ZeroInflatedPoisson, independent of R's pscl.

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import special, stats
from statsmodels.discrete.count_model import ZeroInflatedPoisson
from statsmodels.tools.numdiff import approx_hess3

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any() and (d.visits >= 0).all()

X = sm.add_constant(d[["treated"]].astype(float))
# PINNED: an intercept-only zero-inflation model, passed explicitly. ZeroInflatedPoisson's default
# exog_infl is also intercept-only, the OPPOSITE of pscl's default (same covariates in both parts);
# naming it keeps the two files fitting one model. inflation="logit" named too.
Z = np.ones((len(d), 1))
fit = ZeroInflatedPoisson(d.visits, X, exog_infl=Z, inflation="logit").fit(
    method="newton", maxiter=200, tol=1e-12, disp=False)
assert fit.mle_retvals["converged"], "the zero-inflated fit did not converge"

# PINNED: standard errors from the NUMERICAL Hessian of the model's own log-likelihood, not fit.bse.
# In statsmodels 0.14.2, ZeroInflatedPoisson's bse come from its analytic model.hessian, which at the
# MLE disagrees with the numerical Hessian of the same log-likelihood: on this fixture it gives SEs of
# 0.0956 / 0.0478 / 0.0752 where the log-likelihood's curvature (and R's pscl) gives 0.1060 / 0.0492 /
# 0.0761. fit.bse would understate the uncertainty.
names = list(fit.params.index)
cov = np.linalg.inv(-approx_hess3(fit.params.to_numpy(), fit.model.loglike))
se_all = pd.Series(np.sqrt(np.diag(cov)), index=names)
# statsmodels orders the parameters inflation first, then count.
b, se = fit.params["treated"], se_all["treated"]
infl = fit.params[names[0]]
z = stats.norm.ppf(0.975)
print(fit.params)

print("\n--- HARNESS ---")
print(f"log_rate_ratio={b:.10f}\nlog_rate_ratio_se={se:.10f}")
print(f"log_rate_ratio_lcl={b - z * se:.10f}\nlog_rate_ratio_ucl={b + z * se:.10f}")
print(f"count_intercept={fit.params['const']:.10f}")
print(f"inflation_logit={infl:.10f}")
print(f"structural_zero_prob={special.expit(infl):.10f}")
print(f"n={len(d)}")
