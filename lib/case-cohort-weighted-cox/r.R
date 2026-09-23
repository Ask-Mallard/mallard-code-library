# Case-cohort analysis: a weighted Cox model with a robust variance
#
# In a case-cohort study the exposure is measured in a random subcohort and in every case. The
# subcohort is weighted to stand for the whole cohort (Lin-Ying weights here), and the variance must
# account for the sampling. An ordinary Cox model on the case-cohort sample treats it as the cohort.

library(survival)

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d), !any(duplicated(d$id)), all(d$subcohort == 1 | d$event == 1))

# PINNED: cohort.size (the FULL cohort's size, not the sample's) and method = "LinYing". cch's default
# method is "Prentice"; LinYing uses all subcohort information and its robust variance. The subcohort
# must be a simple random sample of the cohort for these weights.
fit <- cch(Surv(time, event) ~ exposed, data = d, subcoh = ~subcohort, id = ~id,
           cohort.size = 4000, method = "LinYing")

# cch returns an UNNAMED coefficient vector; with one covariate it is position 1.
stopifnot(length(coef(fit)) == 1)
b <- coef(fit)[1]; se <- sqrt(fit$var[1, 1])
z <- qnorm(0.975)
cat(sprintf("hazard ratio %.3f, 95%% %.3f to %.3f\n", exp(b), exp(b - z * se), exp(b + z * se)))

cat("\n--- HARNESS ---\n")
cat(sprintf("log_hazard_ratio=%.10f\nlog_hazard_ratio_se=%.10f\n", b, se))
cat(sprintf("log_hazard_ratio_lcl=%.10f\nlog_hazard_ratio_ucl=%.10f\n", b - z * se, b + z * se))
cat(sprintf("cases=%d\nsubcohort=%d\nsample=%d\n", sum(d$event), sum(d$subcohort), nrow(d)))
