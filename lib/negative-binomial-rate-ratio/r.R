# Negative binomial regression: a rate ratio for overdispersed counts over unequal follow-up
#
# Counts whose variance exceeds their mean (some people prone to events, some not). A Poisson model on
# such data gives the right rate ratio but a standard error that is too small; the negative binomial
# (NB2) model estimates the extra variance through theta. Follow-up enters as an offset, so the ratio
# compares RATES, not counts.

library(MASS)

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d), all(d$followup_years > 0), all(d$events >= 0))

# PINNED: offset(log(followup_years)). Without it the model compares counts, and people followed
# longer look worse.
fit <- glm.nb(events ~ treated + offset(log(followup_years)), data = d,
              control = glm.control(epsilon = 1e-12, maxit = 100))
stopifnot(fit$converged)

est <- coef(summary(fit))
z <- qnorm(0.975)
b <- est["treated", "Estimate"]; se <- est["treated", "Std. Error"]
cat(sprintf("rate ratio %.3f, 95%% %.3f to %.3f; theta %.3f\n", exp(b), exp(b - z * se), exp(b + z * se), fit$theta))

cat("\n--- HARNESS ---\n")
cat(sprintf("log_rate_ratio=%.10f\nlog_rate_ratio_se=%.10f\n", b, se))
cat(sprintf("log_rate_ratio_lcl=%.10f\nlog_rate_ratio_ucl=%.10f\n", b - z * se, b + z * se))
cat(sprintf("theta=%.10f\n", fit$theta))
cat(sprintf("n=%d\n", nrow(d)))
