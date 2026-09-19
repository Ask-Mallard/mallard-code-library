"""Marginal logistic regression by GEE, on binary outcomes clustered within clinic.

The same defaults the R file pins, from the other side:

1. cov_struct=Exchangeable(). statsmodels DEFAULTS TO Independence(), so `sm.GEE(...).fit()` with
   no cov_struct silently fits a different working correlation from the one this file pins.

2. groups= does the clustering, and unlike geepack statsmodels does not need the rows sorted. The
   sort in the R file is not paranoia about this one: it is required there and not here, which is
   itself the kind of difference a translated file loses.

3. Both standard errors, under BOTH working correlations. The gap between the sandwich and the
   model-based error depends on which working correlation produced it, so reporting the pair under
   one structure says nothing about the other -- and the dangerous combination is the one no single
   package defaults to.

standard_errors() returns a bare array rather than a Series, so it is indexed by POSITION taken
from params.index. Indexing it by name raises IndexError, which is how the first run of this file
failed: `.bse` is a Series and `standard_errors()` is not, and the two read identically in code.
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.genmod.cov_struct import Exchangeable, Independence

d = pd.read_csv("fixture.csv")

FORMULA = "outcome ~ exposed + x"


def gee(cov_struct):
    model = sm.GEE.from_formula(FORMULA, groups="clinic", data=d,
                                family=sm.families.Binomial(), cov_struct=cov_struct)
    return model, model.fit()


exch_model, exch = gee(Exchangeable())
ind_model, ind = gee(Independence())

i = list(exch.params.index).index("exposed")
j = list(ind.params.index).index("exposed")

alpha = exch_model.cov_struct.dep_params
alpha = float(np.asarray(alpha).ravel()[0])

print("--- HARNESS ---")
print(f"exposure_log_or={exch.params['exposed']:.10f}")
print(f"exposure_se_robust={np.asarray(exch.standard_errors(cov_type='robust'))[i]:.10f}")
print(f"exposure_se_naive={np.asarray(exch.standard_errors(cov_type='naive'))[i]:.10f}")
print(f"covariate_beta={exch.params['x']:.10f}")
print(f"exposure_log_or_independence={ind.params['exposed']:.10f}")
print(f"exposure_se_robust_independence={np.asarray(ind.standard_errors(cov_type='robust'))[j]:.10f}")
print(f"exposure_se_naive_independence={np.asarray(ind.standard_errors(cov_type='naive'))[j]:.10f}")
print(f"alpha_exchangeable={alpha:.10f}")
