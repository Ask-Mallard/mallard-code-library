# ANCOVA: the treatment effect on a follow-up measurement, adjusted for its baseline value
#
# The Python equivalent of r.R, using statsmodels OLS with its default (model-based) covariance.

import pandas as pd
import statsmodels.formula.api as smf

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any()

# PINNED: follow-up as the outcome, baseline as a covariate; cov_type="nonrobust" named, matching R's
# lm. statsmodels' t-based intervals match R's confint for OLS.
fit = smf.ols("followup ~ treated + baseline", data=d).fit(cov_type="nonrobust")
ci = fit.conf_int(alpha=0.05)

print(fit.summary().tables[1])

print("\n--- HARNESS ---")
print(f"effect={fit.params['treated']:.10f}")
print(f"effect_se={fit.bse['treated']:.10f}")
print(f"effect_lcl={ci.loc['treated', 0]:.10f}\neffect_ucl={ci.loc['treated', 1]:.10f}")
print(f"baseline_coef={fit.params['baseline']:.10f}")
print(f"n={len(d)}")
