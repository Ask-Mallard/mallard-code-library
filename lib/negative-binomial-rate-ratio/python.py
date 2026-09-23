"""Negative binomial regression: a rate ratio for overdispersed counts over unequal follow-up.

Two steps, mirroring what R's glm.nb does: the discrete NegativeBinomial model estimates the
dispersion by maximum likelihood, then a GLM with that dispersion FIXED gives the coefficients and
their standard errors from the expected information, as glm.nb reports them.
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.discrete.discrete_model import NegativeBinomial

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any() and (d.followup_years > 0).all()

X = sm.add_constant(d[["treated"]].astype(float))
offset = np.log(d.followup_years)

# PINNED: the dispersion is ESTIMATED. sm.families.NegativeBinomial() on its own fixes alpha = 1
# (theta = 1) whatever the data, and reports standard errors for that assumed dispersion.
ml = NegativeBinomial(d.events, X, offset=offset, loglike_method="nb2").fit(method="newton", disp=0, maxiter=200, tol=1e-12)
assert ml.mle_retvals["converged"], "the negative binomial fit did not converge"
alpha = float(ml.params["alpha"])
fit = sm.GLM(d.events, X, offset=offset, family=sm.families.NegativeBinomial(alpha=alpha)).fit()

z = stats.norm.ppf(0.975)
b, se = fit.params["treated"], fit.bse["treated"]
print(f"rate ratio {np.exp(b):.3f}, 95% {np.exp(b - z * se):.3f} to {np.exp(b + z * se):.3f}; theta {1 / alpha:.3f}")

print("\n--- HARNESS ---")
print(f"log_rate_ratio={b:.10f}\nlog_rate_ratio_se={se:.10f}")
print(f"log_rate_ratio_lcl={b - z * se:.10f}\nlog_rate_ratio_ucl={b + z * se:.10f}")
print(f"theta={1 / alpha:.10f}")
print(f"n={len(d)}")
