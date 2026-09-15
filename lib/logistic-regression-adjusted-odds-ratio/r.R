# Logistic regression for an ADJUSTED ODDS RATIO, with profile-likelihood and Wald intervals
#
# The coefficient on `exposed` is a log ODDS RATIO, conditional on the covariates in the model.
# It is NOT a risk ratio, and because the odds ratio is non-collapsible it is not the odds ratio
# you would get from the same data with a different covariate set -- even with no confounding.
# See README.md, which prints all three numbers from the generating model.
#
# Base R only. No package is loaded, and none is needed.

d <- read.csv("fixture.csv")

# THE OUTCOME CODING IS CHECKED RATHER THAN ASSUMED. glm(family = binomial) models the probability
# of the SECOND level of a factor, so an outcome arriving as factor(c("yes","no")) is modelled as
# P(yes) only by the accident of alphabetical order -- levels are sorted, "no" sorts first, and the
# odds ratios inverted quietly. A 0/1 numeric column has no such ambiguity, so this file requires
# one and says so instead of trusting the caller.
stopifnot(is.numeric(d$outcome), all(d$outcome %in% c(0, 1)))

# PINNED DEFAULT: link = "logit" is written out. binomial() does default to logit, but the family
# also offers probit, cloglog and log -- and only the logit parameterises an ODDS ratio. A file that
# says `binomial` alone is one edit away from estimating something that is not an odds ratio while
# the variable names and the write-up still say it is.
fit <- glm(outcome ~ exposed + covariate, family = binomial(link = "logit"), data = d)

# THE PINNED DEFAULT THIS ENTRY EXISTS FOR, AND IT IS A DIFFERENCE BETWEEN LANGUAGES RATHER THAN
# INSIDE ONE.
#
#   confint(fit)          PROFILE-LIKELIHOOD intervals. R's method for a glm.
#   confint.default(fit)  WALD intervals, estimate +/- 1.96 SE.
#
# statsmodels' .conf_int() is WALD. So the obvious call in each language -- the one a reader would
# write without thinking about it -- returns a DIFFERENT KIND OF INTERVAL in R than in Python, with
# nothing in either output saying so. They agree only as the two approximations converge, and they
# are furthest apart exactly where it matters: small samples, rare events, near-separation.
#
# Both are reported here, from both languages, and the harness compares both. Neither is "the"
# answer: profile likelihood is the better-behaved default, Wald is what a fitted-model summary
# already prints.
#
# A GUARD, NOT AN ASSUMPTION. confint.glm and profile.glm moved from MASS into stats in R 4.4.0.
# If the method were ever missing, confint() would dispatch to confint.default and return WALD
# intervals under the name `profile` -- silently, in a file whose whole subject is that the two
# differ. That would be this library's own failure mode committed in the entry about it.
if (is.null(getS3method("confint", "glm", optional = TRUE))) {
  stop("no confint.glm method: confint() would silently return Wald intervals labelled profile")
}

ci_profile <- confint(fit, parm = "exposed", level = 0.95)
ci_wald <- confint.default(fit, parm = "exposed", level = 0.95)

cat("adjusted odds ratio for `exposed`, conditional on `covariate`\n")
print(summary(fit)$coefficients)
cat(sprintf("\nodds ratio: %.4f\n", exp(unname(coef(fit)["exposed"]))))
cat(sprintf("  profile-likelihood 95%% CI: %.4f to %.4f\n", exp(ci_profile[1]), exp(ci_profile[2])))
cat(sprintf("  Wald 95%% CI:               %.4f to %.4f\n", exp(ci_wald[1]), exp(ci_wald[2])))
cat("\nThis is an ODDS RATIO. It is not a risk ratio, and it is not collapsible.\n")

cat("\n--- HARNESS ---\n")
cat(sprintf("exposure_log_or=%.10f\n", unname(coef(fit)["exposed"])))
cat(sprintf("exposure_or=%.10f\n", exp(unname(coef(fit)["exposed"]))))
cat(sprintf("exposure_se=%.10f\n", unname(summary(fit)$coefficients["exposed", "Std. Error"])))
cat(sprintf("covariate_beta=%.10f\n", unname(coef(fit)["covariate"])))
cat(sprintf("profile_lo=%.10f\n", unname(ci_profile[1])))
cat(sprintf("profile_hi=%.10f\n", unname(ci_profile[2])))
cat(sprintf("wald_lo=%.10f\n", unname(ci_wald[1])))
cat(sprintf("wald_hi=%.10f\n", unname(ci_wald[2])))
cat(sprintf("n=%d\n", nrow(d)))
