# Nested case-control: conditional logistic regression on risk-set-sampled matched sets.
#
# 1. WHAT THE ODDS RATIO ESTIMATES. When controls are drawn from those still at risk at each case's event
#    time (incidence density sampling), the conditional odds ratio estimates the HAZARD RATIO of the full
#    cohort, with no rare-disease assumption. Here the cohort's hazard ratio is 2.
#
# 2. THE MATCHED SETS MUST STAY TOGETHER. Each set compares a case with the people at risk at that moment;
#    the analysis is conditional on the set (clogit, strata(set)). A person can appear as a control in
#    several sets and later as a case; that is correct and not a duplication to remove.
#
# 3. SAMPLING FROM THE WRONG POOL BIASES IT. Controls must be at risk at the case's time: not everyone
#    who never became a case (cumulative sampling), and not anyone who had already had the event.
#
# 4. clogit's default tie method is an approximation; pinned exact.

library(survival)

d <- read.csv("fixture.csv")
fit <- clogit(case ~ exposed + age + strata(set), data = d, method = "exact")
b <- coef(fit)[["exposed"]]
se <- sqrt(vcov(fit)["exposed", "exposed"])

naive <- coef(glm(case ~ exposed + age, data = d, family = binomial))[["exposed"]]

z <- qnorm(0.975)
cat(sprintf("hazard ratio %.2f (95%% CI %.2f to %.2f) from %d matched sets\n", exp(b), exp(b - z * se),
            exp(b + z * se), length(unique(d$set))))

cat("\n--- HARNESS ---\n")
cat(sprintf("log_hr=%.10f\nlog_hr_se=%.10f\nlog_hr_lcl=%.10f\nlog_hr_ucl=%.10f\n", b, se, b - z * se, b + z * se))
cat(sprintf("age_beta=%.10f\n", coef(fit)[["age"]]))
cat(sprintf("unconditional_log_or=%.10f\nmatched_sets=%d\n", naive, length(unique(d$set))))
