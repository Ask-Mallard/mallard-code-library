# Adjusted risk difference from a binomial model with an identity link
#
# The Python equivalent of r.R, using statsmodels' GLM, an implementation independent of R's glm.

import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

# The identity link is not canonical for the binomial and statsmodels warns so on every fit; that is
# expected here, and the fitted risks are checked to lie inside (0, 1) below.
warnings.filterwarnings("ignore", category=sm.tools.sm_exceptions.DomainWarning)

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any()

# PINNED: link=Identity(); statsmodels' Binomial defaults to the logit link. Start values and a
# tight tolerance as in r.R.
fit = smf.glm("outcome ~ treated + I(age - 60)", data=d,
              family=sm.families.Binomial(link=sm.families.links.Identity())
              ).fit(start_params=np.array([d.outcome.mean(), 0.0, 0.0]), tol=1e-12, maxiter=100)
assert fit.converged and ((fit.fittedvalues > 0) & (fit.fittedvalues < 1)).all()

z = stats.norm.ppf(0.975)
rd, se = fit.params["treated"], fit.bse["treated"]
print(fit.summary().tables[1])

print("\n--- HARNESS ---")
print(f"risk_difference={rd:.10f}\nrisk_difference_se={se:.10f}")
print(f"risk_difference_lcl={rd - z * se:.10f}\nrisk_difference_ucl={rd + z * se:.10f}")
print(f"age_slope={fit.params['I(age - 60)']:.10f}")
print(f"n={len(d)}")
