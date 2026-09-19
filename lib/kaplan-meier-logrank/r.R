# Kaplan-Meier estimation and the log-rank test
#
# The Kaplan-Meier estimator uses partially observed follow-up: a patient censored at ten months
# contributes ten months of "still event-free" and nothing after. The log-rank test compares the
# whole curves rather than survival at any one time, which is why it is not a t-test on a
# percentage.

library(survival)

d <- read.csv("fixture.csv")

# PINNED DEFAULT 1: conf.type = "log-log".
#
# R AND PYTHON DISAGREE HERE AND NEITHER SAYS SO IN ITS OUTPUT:
#   R        survfit()          defaults to "log"
#   Python   lifelines          uses exponential Greenwood, which IS the log-log transform
#
# So R's default band and Python's around the same curve are different bands. log-log is
# pinned because it is the one that cannot leave [0, 1]: the plain log transform produces upper
# limits above 1 in the tail, where a survival probability cannot go.
fit <- survfit(Surv(time, event) ~ exposed, data = d, conf.type = "log-log")

cat("Kaplan-Meier, by arm\n")
print(summary(fit)$table)
cat("\nsurvival at 12 months\n")
s12 <- summary(fit, times = 12)
print(data.frame(strata = s12$strata, surv = s12$surv, lower = s12$lower, upper = s12$upper))

# PINNED DEFAULT 2: rho = 0, which IS the log-rank test.
#
# rho = 1 is the Peto-Peto modification, which weights early differences more heavily. The weighted
# and unweighted tests answer different questions, and a plan that names "the log-rank test" means
# this one.
#
# The log-rank has most power against PROPORTIONAL hazards, which is what this fixture generates.
# It is not a general test that two curves differ: where curves CROSS it can return a p near 1
# while the arms are plainly different, and that is a property of the test, not a bug.
lr <- survdiff(Surv(time, event) ~ exposed, data = d, rho = 0)
cat("\nlog-rank test\n")
print(lr)

tab <- summary(fit)$table
cat("\n--- HARNESS ---\n")
cat(sprintf("median_unexposed=%.10f\n", unname(tab["exposed=0", "median"])))
cat(sprintf("median_exposed=%.10f\n", unname(tab["exposed=1", "median"])))
cat(sprintf("surv_unexposed_12=%.10f\n", s12$surv[1]))
cat(sprintf("surv_exposed_12=%.10f\n", s12$surv[2]))
cat(sprintf("ci_lower_unexposed_12=%.10f\n", s12$lower[1]))
cat(sprintf("ci_upper_unexposed_12=%.10f\n", s12$upper[1]))
cat(sprintf("logrank_chisq=%.10f\n", lr$chisq))
cat(sprintf("n_events=%d\n", sum(d$event)))
cat(sprintf("n=%d\n", nrow(d)))
