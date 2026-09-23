# Linear mixed model for repeated measurements: the treatment-by-time interaction
#
# The Python equivalent of r.R, using statsmodels' MixedLM, independent of R's lme4.

import warnings

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any()

# PINNED: reml=True named (the default), matching lmer's REML = TRUE; visit numeric.
# PINNED: pgtol and factr for L-BFGS. At statsmodels' defaults the patient SD stopped 3.8e-5 short of the
# exact REML value (which a balanced design gives in closed form); lme4 reached it to 1e-8.
# statsmodels warns that pgtol and factr are "not used by MixedLM.fit": they are forwarded to the
# optimizer, and the result moves to the exact REML value, so the warning is silenced deliberately.
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", message="Argument (pgtol|factr) not used")
    fit = smf.mixedlm("y ~ treated * visit", d, groups=d["id"]).fit(reml=True, method=["lbfgs"], pgtol=1e-12, factr=10)
assert fit.converged, "the mixed model did not converge"

b, se = fit.params["treated:visit"], fit.bse["treated:visit"]
z = stats.norm.ppf(0.975)
print(f"difference in slope per visit {b:.3f}, Wald 95% {b - z * se:.3f} to {b + z * se:.3f}")

print("\n--- HARNESS ---")
print(f"interaction={b:.10f}\ninteraction_se={se:.10f}")
print(f"interaction_lcl={b - z * se:.10f}\ninteraction_ucl={b + z * se:.10f}")
print(f"time_slope={fit.params['visit']:.10f}")
print(f"sd_patient={float(np.sqrt(fit.cov_re.iloc[0, 0])):.10f}\nsd_residual={float(np.sqrt(fit.scale)):.10f}")
print(f"n={len(d)}\npatients={d.id.nunique()}")
