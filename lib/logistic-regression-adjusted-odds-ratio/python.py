# Logistic regression for an ADJUSTED ODDS RATIO, with profile-likelihood and Wald intervals
#
# The Python counterpart of the R file beside this one. Same model, same two intervals -- but the
# profile interval has to be COMPUTED here, because statsmodels does not provide one for GLM. That
# asymmetry is the entry's subject, not an inconvenience: see the comment above the profile block.

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import optimize, stats

d = pd.read_csv("fixture.csv")

# Checked, not assumed -- the mirror of the R file's guard. statsmodels will happily fit a Binomial
# GLM to a column of "yes"/"no" strings after a cast, or to a 1/2 coding, and the coefficient then
# refers to whichever level the cast made 1.
assert d["outcome"].isin([0, 1]).all(), "outcome must be coded 0/1"

# PINNED DEFAULT: link=Logit() is written out, matching the R file. Binomial() does default to the
# logit, but it also offers probit and cloglog, and only the logit parameterises an ODDS ratio.
model = smf.glm("outcome ~ exposed + covariate",
                data=d,
                family=sm.families.Binomial(link=sm.families.links.Logit()))

# PINNED DEFAULT: the model-based (nonrobust) covariance is CORRECT here and is deliberately left
# alone. This is the opposite of the modified-Poisson entry beside it, where the Poisson variance
# assumption is false by construction and HC0 is mandatory. A correctly specified logistic
# likelihood has the information matrix as its variance; reaching for cov_type="HC0" by habit would
# discard the efficient standard error for no reason and break agreement with the R file.
fit = model.fit()

# WALD, and statsmodels' .conf_int() IS Wald -- estimate +/- z * SE. R's confint() on a glm is NOT:
# it is profile likelihood. So the obvious one-liner in each language returns a different KIND of
# interval, with nothing in either output saying so. That is what this entry pins.
wald = fit.conf_int(alpha=0.05).loc["exposed"]

# PROFILE LIKELIHOOD, COMPUTED EXPLICITLY, because statsmodels has no equivalent of R's
# confint.glm for a GLM. The bounds are the values of the exposure coefficient where the deviance
# rises above its minimum by qchisq(0.95, 1) -- the same criterion R's method uses.
#
# Holding the coefficient at `c` is done with an OFFSET: c * exposed enters the linear predictor
# with its coefficient fixed at 1, and the remaining parameters are re-fitted. That IS the profile
# deviance, not an approximation of it.
#
# WHERE THIS DIFFERS FROM R, AND IT IS THE R SIDE THAT APPROXIMATES: R's confint.glm evaluates the
# profile on a grid and interpolates the crossing with a spline, while brentq solves for the root
# directly. The two therefore disagree by interpolation error rather than by anything statistical.
# MEASURED, because the point of this library is to read the number rather than assert it: the
# largest gap between the two engines on any profile bound is 5.9e-7, so no loosened tolerance is
# needed and none is declared. Swapping a Wald bound in for a profile one moves the same number by
# 2.3e-4, which the shared 1e-4 tolerance still catches.
CRITICAL = stats.chi2.ppf(0.95, df=1)
deviance_hat = fit.deviance
exposed = d["exposed"].to_numpy(dtype=float)


def profile_gap(c: float) -> float:
    """Deviance excess over the critical value with the exposure coefficient held at `c`."""
    restricted = smf.glm("outcome ~ covariate",
                         data=d,
                         family=sm.families.Binomial(link=sm.families.links.Logit()),
                         offset=c * exposed).fit()
    return restricted.deviance - deviance_hat - CRITICAL


estimate = fit.params["exposed"]
se = fit.bse["exposed"]
# The bracket is deliberately wide: the profile bound sits near the Wald one, so +/- 10 SE contains
# it with room to spare, and brentq requires a sign change rather than a good guess.
profile_lo = optimize.brentq(profile_gap, estimate - 10 * se, estimate, xtol=1e-12, rtol=8.9e-16)
profile_hi = optimize.brentq(profile_gap, estimate, estimate + 10 * se, xtol=1e-12, rtol=8.9e-16)

print("adjusted odds ratio for `exposed`, conditional on `covariate`")
print(fit.summary())
print(f"\nodds ratio: {np.exp(estimate):.4f}")
print(f"  profile-likelihood 95% CI: {np.exp(profile_lo):.4f} to {np.exp(profile_hi):.4f}")
print(f"  Wald 95% CI:               {np.exp(wald[0]):.4f} to {np.exp(wald[1]):.4f}")
print("\nThis is an ODDS RATIO. It is not a risk ratio, and it is not collapsible.")

print("\n--- HARNESS ---")
print(f"exposure_log_or={estimate:.10f}")
print(f"exposure_or={np.exp(estimate):.10f}")
print(f"exposure_se={se:.10f}")
print(f"covariate_beta={fit.params['covariate']:.10f}")
print(f"profile_lo={profile_lo:.10f}")
print(f"profile_hi={profile_hi:.10f}")
print(f"wald_lo={wald[0]:.10f}")
print(f"wald_hi={wald[1]:.10f}")
print(f"n={len(d)}")
