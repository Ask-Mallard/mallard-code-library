"""Stepped-wedge trial: the Hussey-Hughes linear mixed model.

The Python equivalent of r.R, fitted with statsmodels' MixedLM (REML). The intervention's standard error
is computed as lme4 reports it: (X' V^-1 X)^-1 with V built from the REML variance components, the
covariance of the fixed effects CONDITIONAL on those components. statsmodels' own bse inverts the joint
Hessian of fixed effects and variance parameters, which in an unbalanced design such as a stepped wedge
gives a slightly different number (0.66797 against 0.66795 here). The conditional form is the standard
one, and the one small-sample corrections such as Kenward-Roger start from.
"""

import warnings

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any()

# PINNED: C(period) as fixed effects, reml=True, and tight L-BFGS settings (statsmodels' defaults stop
# short of the REML optimum; see linear-mixed-model-repeated-measures). The "not used" warning for
# pgtol and factr is silenced because they do reach the optimizer. statsmodels also warns "Random
# effects covariance is singular" when an optimizer ITERATE touches the boundary; that warning is
# silenced too, and the FINAL ward variance is asserted positive below, so a truly singular fit still stops.
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", message="Argument (pgtol|factr) not used")
    warnings.filterwarnings("ignore", message="Random effects covariance is singular")
    fit = smf.mixedlm("y ~ intervention + C(period)", d, groups=d["ward"]).fit(
        reml=True, method=["lbfgs"], pgtol=1e-12, factr=10)
assert fit.converged, "the mixed model did not converge"
assert float(fit.cov_re.iloc[0, 0]) > 1e-8, "the ward variance is at zero: the fit is singular"

# PINNED: the fixed-effect covariance conditional on the REML variance components (see the docstring).
X = fit.model.exog
s2_resid, s2_ward = float(fit.scale), float(fit.cov_re.iloc[0, 0])
information = np.zeros((X.shape[1], X.shape[1]))
for idx in d.groupby("ward").indices.values():
    Xc = X[idx]
    V = s2_resid * np.eye(len(idx)) + s2_ward * np.ones((len(idx), len(idx)))
    information += Xc.T @ np.linalg.solve(V, Xc)
j = list(fit.model.exog_names).index("intervention")
b = float(fit.fe_params["intervention"])
se = float(np.sqrt(np.linalg.inv(information)[j, j]))

z = stats.norm.ppf(0.975)
print(f"intervention effect {b:.3f}, Wald 95% {b - z * se:.3f} to {b + z * se:.3f}")

print("\n--- HARNESS ---")
print(f"effect={b:.10f}\neffect_se={se:.10f}\neffect_lcl={b - z * se:.10f}\neffect_ucl={b + z * se:.10f}")
print(f"sd_ward={np.sqrt(s2_ward):.10f}\nsd_residual={np.sqrt(s2_resid):.10f}")
print(f"n={len(d)}\nwards={d.ward.nunique()}")
