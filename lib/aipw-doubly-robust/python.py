"""Augmented inverse-probability weighting (AIPW): a doubly robust average treatment effect.

The same estimator as r.R: logistic propensity and outcome models, each patient's predicted risks
corrected by the weighted residual of the arm they were in, unstabilized and untrimmed, with a standard
error from the influence function. Inverse-probability weighting alone and g-computation alone are
reported beside it. See r.R for the formula and for what the influence-function SE assumes.
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

d = pd.read_csv("fixture.csv")
a = d["treated"].to_numpy(float)
y = d["event"].to_numpy(float)

ps = smf.glm("treated ~ severity + comorbid", d, family=sm.families.Binomial()).fit(tol=1e-12)
e = np.asarray(ps.fittedvalues)
om = smf.glm("event ~ treated + severity + comorbid", d, family=sm.families.Binomial()).fit(tol=1e-12)
m1 = np.asarray(om.predict(d.assign(treated=1)))
m0 = np.asarray(om.predict(d.assign(treated=0)))

# PINNED: the AIPW estimating equation, unstabilized, untrimmed.
phi1 = m1 + a * (y - m1) / e
phi0 = m0 + (1 - a) * (y - m0) / (1 - e)
ate = phi1.mean() - phi0.mean()
n = len(d)
IF = phi1 - phi0 - ate
se = float(np.sqrt(np.sum(IF ** 2)) / n)

ipw = np.mean(a * y / e) - np.mean((1 - a) * y / (1 - e))
gcomp = m1.mean() - m0.mean()

z = stats.norm.ppf(0.975)
print(f"AIPW risk difference {ate:.4f} (95% CI {ate - z * se:.4f} to {ate + z * se:.4f})")

print("\n--- HARNESS ---")
print(f"ate={ate:.10f}\nate_se={se:.10f}\nate_lcl={ate - z * se:.10f}\nate_ucl={ate + z * se:.10f}")
print(f"ipw_only={ipw:.10f}\ngcomp_only={gcomp:.10f}")
print(f"n={n}")
