"""Self-controlled case series (SCCS): conditional Poisson regression within each case.

The same analysis as r.R: a Poisson regression with one fixed effect per case and log(interval length)
as the offset, which has the same estimate and standard error as the SCCS conditional likelihood. Age
band is in the model; the fit without it is reported beside it. The GLM runs to a tight tolerance so the
400 child effects do not leave the risk coefficient short of R's. See r.R for the design's assumptions.
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

d = pd.read_csv("fixture.csv")
off = np.log(d["days"])
# PINNED: wls_method="qr". The IRLS default solves each step by SVD, which failed to converge ("SVD did
# not converge") on one of the 100 calibration fixtures; QR fits all of them to the same estimates.
fit = smf.glm("events ~ risk + C(age_band) + C(child)", d, family=sm.families.Poisson(),
              offset=off).fit(tol=1e-12, wls_method="qr")
b, se = float(fit.params["risk"]), float(fit.bse["risk"])

no_age = smf.glm("events ~ risk + C(child)", d, family=sm.families.Poisson(),
                 offset=off).fit(tol=1e-12, wls_method="qr").params["risk"]

z = stats.norm.ppf(0.975)
print(f"incidence rate ratio {np.exp(b):.2f} (95% CI {np.exp(b - z * se):.2f} to {np.exp(b + z * se):.2f})")

print("\n--- HARNESS ---")
print(f"log_irr={b:.10f}\nlog_irr_se={se:.10f}\nlog_irr_lcl={b - z * se:.10f}\nlog_irr_ucl={b + z * se:.10f}")
print(f"log_irr_age_ignored={no_age:.10f}")
print(f"cases={d.child.nunique()}\nevents={int(d.events.sum())}")
