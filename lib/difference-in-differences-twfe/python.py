"""Difference-in-differences: two-way fixed effects for a policy adopted at one time.

The same analysis as r.R: hospital and quarter fixed effects, the coefficient on the policy indicator,
and a standard error clustered on hospital. statsmodels' cluster covariance applies the same small-sample
factor as sandwich::vcovCL(type = "HC1"), G/(G-1) x (N-1)/(N-K). See r.R for parallel trends and for why
this estimator is for a single adoption date only.
"""

import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

d = pd.read_csv("fixture.csv")
fit = smf.ols("los ~ policy + C(hospital) + C(quarter)", d).fit(
    cov_type="cluster", cov_kwds={"groups": d["hospital"]})
did, se = float(fit.params["policy"]), float(fit.bse["policy"])

g = lambda a, p: d.loc[(d.adopter == a) & (d.post == p), "los"].mean()
pre_post_adopters = g(1, 1) - g(1, 0)
post_between_groups = g(1, 1) - g(0, 1)

z = stats.norm.ppf(0.975)
print(f"DiD {did:.3f} (95% CI {did - z * se:.3f} to {did + z * se:.3f})")

print("\n--- HARNESS ---")
print(f"did={did:.10f}\ndid_se={se:.10f}\ndid_lcl={did - z * se:.10f}\ndid_ucl={did + z * se:.10f}")
print(f"pre_post_adopters={pre_post_adopters:.10f}\npost_between_groups={post_between_groups:.10f}")
print(f"hospitals={d.hospital.nunique()}\nn={len(d)}")
