# Rate regression on overdispersed counts with unequal follow-up
#
# The Python equivalent of the R file beside this one.

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.discrete.discrete_model import NegativeBinomial

d = pd.read_csv("fixture.csv")

X = sm.add_constant(d[["exposed", "covariate"]])

# PINNED DEFAULT 1: THE OFFSET IS log(person_time), passed as an offset rather than a column.
# statsmodels takes it already logged, which is the same convention as R's offset(log(t)). The
# offset must always be the LOG of the person-time; passing the person-time itself fits a different
# model that still runs.
# Passing person_time where log(person_time) belongs runs, converges, and answers a different
# question.
offset = np.log(d["person_time"].to_numpy())

model = sm.GLM(d["events"], X, family=sm.families.Poisson(), offset=offset)

# PINNED DEFAULT 2: cov_type="HC0", the robust sandwich standard errors.
# A Poisson model assumes variance equals mean; these counts have a variance 2.4 times their mean.
# The coefficients stay consistent under that and the model-based standard errors do not.
res = model.fit(cov_type="HC0")
naive = model.fit()

print("Poisson rate regression, robust standard errors")
print(res.summary())
print("\nmodel-based standard errors, FOR COMPARISON ONLY -- they assume variance equals mean")
print(naive.bse)

# The negative binomial alternative.
#
# STATSMODELS REPORTS alpha AND R REPORTS theta, AND THEY ARE RECIPROCALS. nb2 parameterises the
# variance as mu + alpha*mu^2; R's glm.nb uses mu + mu^2/theta. The harness key below is alpha,
# and the R file converts its theta to alpha rather than the two files reporting different
# quantities under one name.
nb = NegativeBinomial(d["events"], X, loglike_method="nb2",
                      offset=offset).fit(method="bfgs", maxiter=1000, disp=0)
print("\nnegative binomial")
print(nb.summary())

print("\n--- HARNESS ---")
print(f"exposure_log_rr={res.params['exposed']:.10f}")
print(f"exposure_se={res.bse['exposed']:.10f}")
print(f"exposure_se_naive={naive.bse['exposed']:.10f}")
print(f"covariate_beta={res.params['covariate']:.10f}")
print(f"nb_exposure_log_rr={nb.params['exposed']:.10f}")
print(f"nb_alpha={nb.params['alpha']:.10f}")
print(f"events={int(d['events'].sum())}")
print(f"n={len(d)}")
