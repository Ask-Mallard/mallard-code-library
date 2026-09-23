# Cause-specific Cox regression under competing risks
#
# The hazard of relapse among people still alive and relapse-free, with deaths treated as censoring
# FOR THIS MODEL ONLY. That is correct for the cause-specific hazard (a rate), and it is why this
# hazard ratio does not translate into the proportion who relapse: that also depends on how fast the
# competing death removes people. Report the cumulative incidence (competing-risks-cumulative-incidence)
# beside it when the question is "how many relapse".

library(survival)

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d), all(d$status %in% c(0, 1, 2)))

# PINNED: the event is status == 1 only. Counting deaths as events (a composite) answers a different
# question, and dropping the people who died removes person-time that was at risk of relapse.
fit <- coxph(Surv(time, status == 1) ~ treated, data = d, ties = "efron")

z <- qnorm(0.975)
b <- coef(fit)[["treated"]]; se <- sqrt(vcov(fit)["treated", "treated"])
cat(sprintf("cause-specific HR for relapse %.3f, 95%% %.3f to %.3f\n", exp(b), exp(b - z * se), exp(b + z * se)))

cat("\n--- HARNESS ---\n")
cat(sprintf("log_cs_hazard_ratio=%.10f\nlog_cs_hazard_ratio_se=%.10f\n", b, se))
cat(sprintf("log_cs_hazard_ratio_lcl=%.10f\nlog_cs_hazard_ratio_ucl=%.10f\n", b - z * se, b + z * se))
cat(sprintf("relapses=%d\ncompeting_deaths=%d\nn=%d\n", sum(d$status == 1), sum(d$status == 2), nrow(d)))
