# Multinomial logistic regression for an unordered categorical outcome
#
# The Python equivalent of r.R, using statsmodels' MNLogit, independent of R's nnet::multinom.
# MNLogit takes the LOWEST integer code as the reference, so the codes are assigned explicitly.

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any()

# PINNED: home = 0 is the reference; the codes are written out, not taken from sorted labels.
levels = ["home", "rehab", "nursing"]
y = d.destination.map({k: i for i, k in enumerate(levels)})
assert not y.isna().any()

X = sm.add_constant(d[["treated"]].astype(float))
fit = sm.MNLogit(y, X).fit(method="newton", maxiter=200, tol=1e-12, disp=False)
assert fit.mle_retvals["converged"], "the multinomial fit did not converge"

z = stats.norm.ppf(0.975)
print(fit.params)

print("\n--- HARNESS ---")
# Column j of params/bse is the equation for code j + 1 against the reference.
for j, k in enumerate(levels[1:]):
    b, se = fit.params.iloc[1, j], fit.bse.iloc[1, j]
    print(f"{k}_intercept={fit.params.iloc[0, j]:.10f}")
    print(f"{k}_log_rrr={b:.10f}\n{k}_log_rrr_se={se:.10f}")
    print(f"{k}_log_rrr_lcl={b - z * se:.10f}\n{k}_log_rrr_ucl={b + z * se:.10f}")
print(f"n={len(d)}")
