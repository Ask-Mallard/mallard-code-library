# Andersen-Gill model for recurrent events, with a robust (clustered) variance
#
# Uses every admission, not just the first. The rate ratio is population-averaged, and the robust
# variance, clustered on the person, accounts for some people being admitted far more often than
# others. Without it the standard error is too small.

library(survival)

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d), all(d$stop > d$start))

# PINNED: cluster = id, which makes vcov() the robust (Lin-Wei) sandwich. Without it coxph reports the
# model-based variance, which treats every admission as independent.
fit <- coxph(Surv(start, stop, event) ~ treated, data = d, cluster = id, ties = "efron")

b <- coef(fit)[["treated"]]
se <- sqrt(vcov(fit)["treated", "treated"])
naive_se <- sqrt(fit$naive.var[1, 1])
z <- qnorm(0.975)
cat(sprintf("rate ratio %.3f, robust 95%% %.3f to %.3f (robust SE %.4f, naive %.4f)\n",
            exp(b), exp(b - z * se), exp(b + z * se), se, naive_se))

cat("\n--- HARNESS ---\n")
cat(sprintf("log_rate_ratio=%.10f\nlog_rate_ratio_robust_se=%.10f\n", b, se))
cat(sprintf("log_rate_ratio_lcl=%.10f\nlog_rate_ratio_ucl=%.10f\n", b - z * se, b + z * se))
cat(sprintf("log_rate_ratio_naive_se=%.10f\n", naive_se))
cat(sprintf("events=%d\npeople=%d\n", sum(d$event), length(unique(d$id))))
