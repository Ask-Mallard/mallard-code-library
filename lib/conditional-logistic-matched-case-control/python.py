# Conditional logistic regression for a 1:m matched case-control study
#
# The Python equivalent of survival::clogit. statsmodels calls it ConditionalLogit, and it
# conditions on the group rather than estimating a coefficient per stratum.
#
# When the matched sets are risk sets (controls drawn from those still at risk at each case's event
# time, as in a nested case-control study), this conditional odds ratio estimates the cohort's hazard
# (rate) ratio, with no rare-disease assumption.

import numpy as np
import pandas as pd
from statsmodels.discrete.conditional_models import ConditionalLogit

d = pd.read_csv("fixture.csv")

# PINNED DEFAULT: no intercept. A conditional likelihood has no intercept to estimate, because the
# stratum terms are conditioned away. statsmodels does not add one here, unlike the formula API
# used elsewhere, so the design matrix is built explicitly rather than from a formula string.
exog = d[["exposed", "covariate"]].astype(float)
endog = d["case"].astype(int)
groups = d["set_id"].astype(int)

model = ConditionalLogit(endog, exog, groups=groups)
result = model.fit(disp=0)

print("estimates")
print(result.summary())

print("\n--- HARNESS ---")
print(f"exposure_log_or={result.params['exposed']:.10f}")
print(f"exposure_se={result.bse['exposed']:.10f}")
print(f"covariate_beta={result.params['covariate']:.10f}")
print(f"n_sets={groups.nunique()}")
