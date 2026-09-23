# Weibull regression: an accelerated failure time (time ratio), and the equivalent hazard ratio
#
# A parametric survival model assumes a shape for the hazard (Weibull: rising or falling steadily). In
# return it gives the time ratio ("treatment stretches survival times by 49%") and allows
# extrapolation beyond follow-up, which Cox regression does not; extrapolation is only as good as the
# assumed shape.

library(survival)

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d))

fit <- survreg(Surv(time, event) ~ treated, data = d, dist = "weibull")

# PINNED: survreg's coefficients are on the log TIME scale (AFT), and its "scale" is 1 / shape. A
# positive coefficient means LONGER survival. The Weibull hazard ratio is exp(-coefficient / scale);
# reading the coefficient as a log hazard ratio gets the sign and the size wrong.
b <- coef(fit)[["treated"]]; se <- sqrt(vcov(fit)["treated", "treated"])
shape <- 1 / fit$scale
z <- qnorm(0.975)
cat(sprintf("time ratio %.3f, 95%% %.3f to %.3f; shape %.3f; hazard ratio %.3f\n",
            exp(b), exp(b - z * se), exp(b + z * se), shape, exp(-b * shape)))

cat("\n--- HARNESS ---\n")
cat(sprintf("log_time_ratio=%.10f\nlog_time_ratio_se=%.10f\n", b, se))
cat(sprintf("log_time_ratio_lcl=%.10f\nlog_time_ratio_ucl=%.10f\n", b - z * se, b + z * se))
cat(sprintf("shape=%.10f\n", shape))
cat(sprintf("log_hazard_ratio=%.10f\n", -b * shape))
cat(sprintf("n=%d\n", nrow(d)))
