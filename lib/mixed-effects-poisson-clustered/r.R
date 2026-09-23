# Mixed-effects Poisson regression for clustered counts
#
# Infection counts per patient with length of stay as the exposure, patients clustered in wards. A random
# intercept per ward carries the between-ward variation; ignoring it gives a standard error that is too
# small. With a log link the rate ratio is collapsible over the random intercept, so the ward-conditional
# rate ratio reported here is also the population-average one.

suppressPackageStartupMessages(library(lme4))

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d), all(d$days > 0), all(d$infections >= 0))

# PINNED: offset(log(days)), so the ratio compares RATES per patient-day; family = poisson; and the
# Laplace approximation (nAGQ = 1), lme4's default, named.
fit <- glmer(infections ~ intervention + offset(log(days)) + (1 | ward), data = d, family = poisson, nAGQ = 1)

co <- summary(fit)$coefficients
b <- co["intervention", "Estimate"]; se <- co["intervention", "Std. Error"]
z <- qnorm(0.975)
sd_ward <- attr(VarCorr(fit)$ward, "stddev")[[1]]
cat(sprintf("rate ratio %.3f, Wald 95%% %.3f to %.3f; ward SD %.3f\n", exp(b), exp(b - z * se), exp(b + z * se), sd_ward))

cat("\n--- HARNESS ---\n")
cat(sprintf("log_rate_ratio=%.10f\nlog_rate_ratio_se=%.10f\n", b, se))
cat(sprintf("log_rate_ratio_lcl=%.10f\nlog_rate_ratio_ucl=%.10f\n", b - z * se, b + z * se))
cat(sprintf("sd_ward=%.10f\n", sd_ward))
cat(sprintf("wards=%d\nn=%d\n", length(unique(d$ward)), nrow(d)))
