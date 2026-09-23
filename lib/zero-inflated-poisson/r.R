# Zero-inflated Poisson regression for counts with structural zeros
#
# Two parts: a logistic model for being a structural zero (never at risk) and a Poisson model for the
# count among those at risk. The count part's treatment coefficient is a log rate ratio AMONG THOSE AT
# RISK, not in the whole population; the report must say which. A zero-inflated model is justified
# by a mechanism that produces structural zeros, not by "there are many zeros": a negative binomial
# often fits excess zeros too (see negative-binomial-rate-ratio).

suppressPackageStartupMessages(library(pscl))

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d), all(d$visits >= 0))

# PINNED: the formula's second part, "| 1", gives the zero-inflation model its own (intercept-only)
# covariates. Without the bar, zeroinfl uses the SAME covariates in both parts. PINNED: dist =
# "poisson" (the default) and a tight optimizer tolerance.
fit <- zeroinfl(visits ~ treated | 1, data = d, dist = "poisson",
                control = zeroinfl.control(reltol = 1e-14, maxit = 10000))
stopifnot(fit$converged)

s <- summary(fit)$coefficients
z <- qnorm(0.975)
b <- s$count["treated", "Estimate"]; se <- s$count["treated", "Std. Error"]
infl <- s$zero["(Intercept)", "Estimate"]
print(s)

cat("\n--- HARNESS ---\n")
cat(sprintf("log_rate_ratio=%.10f\nlog_rate_ratio_se=%.10f\n", b, se))
cat(sprintf("log_rate_ratio_lcl=%.10f\nlog_rate_ratio_ucl=%.10f\n", b - z * se, b + z * se))
cat(sprintf("count_intercept=%.10f\n", s$count["(Intercept)", "Estimate"]))
cat(sprintf("inflation_logit=%.10f\n", infl))
cat(sprintf("structural_zero_prob=%.10f\n", plogis(infl)))
cat(sprintf("n=%d\n", nrow(d)))
