# Cox regression with a time-varying exposure, in counting-process (start, stop] form
#
# When exposure begins during follow-up, each person's time is split at the start, and the exposure
# is whatever it was during each interval. The hazard ratio then compares CURRENTLY exposed with
# currently unexposed person-time. Coding "ever exposed" at baseline instead creates immortal time
# and biases the estimate toward benefit.

library(survival)

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d), all(d$stop > d$start))

# PINNED: Surv(start, stop, event), the counting-process form, and ties = "efron" named (it is
# coxph's default; lifelines uses Efron too).
fit <- coxph(Surv(start, stop, event) ~ exposed, data = d, ties = "efron")

z <- qnorm(0.975)
b <- coef(fit)[["exposed"]]; se <- sqrt(vcov(fit)["exposed", "exposed"])
cat(sprintf("hazard ratio for current exposure %.3f, 95%% %.3f to %.3f\n", exp(b), exp(b - z * se), exp(b + z * se)))

cat("\n--- HARNESS ---\n")
cat(sprintf("log_hazard_ratio=%.10f\nlog_hazard_ratio_se=%.10f\n", b, se))
cat(sprintf("log_hazard_ratio_lcl=%.10f\nlog_hazard_ratio_ucl=%.10f\n", b - z * se, b + z * se))
cat(sprintf("events=%d\npeople=%d\n", sum(d$event), length(unique(d$id))))
