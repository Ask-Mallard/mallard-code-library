"""Nested case-control: conditional logistic regression on risk-set-sampled matched sets.

The same analysis as r.R with statsmodels' ConditionalLogit, grouped on the matched set, no intercept,
and method="newton" (its default optimizer stops short of the maximum; see
case-crossover-conditional-logistic). The unconditional logistic regression that ignores the sets is
reported beside it. See r.R for why the conditional odds ratio estimates the cohort's hazard ratio.
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.discrete.conditional_models import ConditionalLogit

d = pd.read_csv("fixture.csv")
fit = ConditionalLogit(d["case"].astype(int), d[["exposed", "age"]].astype(float),
                       groups=d["set"].astype(int)).fit(disp=0, method="newton", tol=1e-12)
b, se = float(fit.params["exposed"]), float(fit.bse["exposed"])

naive = sm.Logit(d["case"], sm.add_constant(d[["exposed", "age"]].astype(float))).fit(disp=0).params["exposed"]

z = stats.norm.ppf(0.975)
print(f"hazard ratio {np.exp(b):.2f} (95% CI {np.exp(b - z * se):.2f} to {np.exp(b + z * se):.2f}) "
      f"from {d['set'].nunique()} matched sets")

print("\n--- HARNESS ---")
print(f"log_hr={b:.10f}\nlog_hr_se={se:.10f}\nlog_hr_lcl={b - z * se:.10f}\nlog_hr_ucl={b + z * se:.10f}")
print(f"age_beta={fit.params['age']:.10f}")
print(f"unconditional_log_or={naive:.10f}\nmatched_sets={d['set'].nunique()}")
