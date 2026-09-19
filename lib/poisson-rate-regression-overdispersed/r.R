# Rate regression on overdispersed counts with unequal follow-up
#
# Poisson regression with an offset estimates a RATE ratio: events per person-time, not events per
# patient. With follow-up ranging from six months to five years those are different questions, and
# only the offset makes it the first one.

library(sandwich)
library(lmtest)
library(MASS)

d <- read.csv("fixture.csv")

# PINNED DEFAULT 1: THE OFFSET IS log(person_time), AND IT IS A FIXED TERM, NOT A COVARIATE.
# offset() forces the coefficient to exactly 1 rather than estimating it. Writing
# `+ log(person_time)` instead fits a model that runs, converges, and answers a different question.
fit <- glm(events ~ exposed + covariate + offset(log(person_time)),
           family = poisson, data = d)

# PINNED DEFAULT 2: THE STANDARD ERRORS ARE ROBUST, NOT THE MODEL-BASED ONES.
# A Poisson model assumes variance equals mean. These counts have a variance 2.4 times their mean,
# and under that the COEFFICIENTS ARE STILL CONSISTENT while the model-based standard errors are
# materially too small. The estimate is right and the interval is too narrow, which is the failure
# that matters: an analysis that is confidently wrong rather than visibly wrong.
V <- vcovHC(fit, type = "HC0")

cat("Poisson rate regression, robust standard errors\n")
print(coeftest(fit, vcov. = V))
cat("\nmodel-based standard errors, FOR COMPARISON ONLY -- they assume variance equals mean\n")
print(summary(fit)$coefficients)
cat(sprintf("\nresidual deviance %.2f on %d df -- far above its df, which is the overdispersion\n",
            fit$deviance, fit$df.residual))

# The negative binomial alternative: it models the overdispersion rather than correcting the
# standard errors for it, and gives a fully specified likelihood.
#
# R REPORTS theta AND statsmodels REPORTS alpha, AND THEY ARE RECIPROCALS. glm.nb
# parameterises the variance as mu + mu^2/theta; statsmodels uses mu + alpha*mu^2. Reading R's
# theta of 2 as if it were alpha would describe data four times less dispersed than it is. The
# fixture's gamma variance is 0.5, so theta is 2 and alpha is 0.5.
nb <- glm.nb(events ~ exposed + covariate + offset(log(person_time)), data = d)
cat("\nnegative binomial\n")
print(summary(nb)$coefficients)
cat(sprintf("theta = %.6f  (the same thing other packages call alpha = %.6f)\n",
            nb$theta, 1 / nb$theta))

cat("\n--- HARNESS ---\n")
cat(sprintf("exposure_log_rr=%.10f\n", unname(coef(fit)["exposed"])))
cat(sprintf("exposure_se=%.10f\n", sqrt(V["exposed", "exposed"])))
cat(sprintf("exposure_se_naive=%.10f\n", summary(fit)$coefficients["exposed", "Std. Error"]))
cat(sprintf("covariate_beta=%.10f\n", unname(coef(fit)["covariate"])))
cat(sprintf("nb_exposure_log_rr=%.10f\n", unname(coef(nb)["exposed"])))
cat(sprintf("nb_alpha=%.10f\n", 1 / nb$theta))
cat(sprintf("events=%d\n", sum(d$events)))
cat(sprintf("n=%d\n", nrow(d)))
