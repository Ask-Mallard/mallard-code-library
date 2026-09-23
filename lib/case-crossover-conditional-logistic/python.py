"""Case-crossover: each patient's hazard window against their own referent windows.

The same analysis as r.R with statsmodels' ConditionalLogit, which conditions on the patient (the group)
rather than estimating an intercept per patient. With one case window per stratum the conditional
likelihood is exact, which is what clogit(method = "exact") computes. The ordinary logistic regression
that ignores the matching is reported beside it. See r.R for the design's assumptions.
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.discrete.conditional_models import ConditionalLogit

d = pd.read_csv("fixture.csv")

# No intercept: the patient terms are conditioned away.
# PINNED: method="newton". ConditionalLogit's default optimizer (BFGS) stopped 5.4e-5 short of the
# maximum on this fixture, with gtol=1e-10 too; Newton reaches the exact conditional MLE that clogit does.
fit = ConditionalLogit(d["hazard_window"].astype(int), d[["exposed"]].astype(float),
                       groups=d["patient"].astype(int)).fit(disp=0, method="newton", tol=1e-12)
b, se = float(fit.params["exposed"]), float(fit.bse["exposed"])

naive = sm.Logit(d["hazard_window"], sm.add_constant(d[["exposed"]].astype(float))).fit(disp=0).params["exposed"]
informative = int((d.groupby("patient")["exposed"].nunique() > 1).sum())

z = stats.norm.ppf(0.975)
print(f"odds ratio {np.exp(b):.2f} (95% CI {np.exp(b - z * se):.2f} to {np.exp(b + z * se):.2f}); "
      f"{informative} of {d.patient.nunique()} patients informative")

print("\n--- HARNESS ---")
print(f"log_or={b:.10f}\nlog_or_se={se:.10f}\nlog_or_lcl={b - z * se:.10f}\nlog_or_ucl={b + z * se:.10f}")
print(f"unconditional_log_or={naive:.10f}\ninformative_patients={informative}\nn_patients={d.patient.nunique()}")
