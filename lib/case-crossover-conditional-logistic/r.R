# Case-crossover: each patient's hazard window against their own referent windows.
#
# 1. EVERY PATIENT IS THEIR OWN CONTROL. Only cases are sampled; exposure just before the event is
#    compared with exposure in referent windows of the same person. Anything that does not change over
#    those windows (sex, genes, usual habits) is removed by design, measured or not. It suits transient
#    exposures with acute effects.
#
# 2. THE ANALYSIS IS CONDITIONAL LOGISTIC REGRESSION, STRATIFIED ON THE PATIENT. An ordinary logistic
#    regression ignores the matching and is biased; here it gives a smaller odds ratio.
#
# 3. REFERENT WINDOWS: TIME-STRATIFIED. Referents in fixed strata (for example the same weekday in the
#    same month) avoid the overlap bias of choosing referents relative to the event date. Time trends in
#    exposure within a stratum still bias the estimate.
#
# 4. clogit's default tie method is Breslow's approximation passed through coxph; pinned exact, as in
#    conditional-logistic-matched-case-control.

library(survival)

d <- read.csv("fixture.csv")
fit <- clogit(hazard_window ~ exposed + strata(patient), data = d, method = "exact")
b <- coef(fit)[["exposed"]]
se <- sqrt(vcov(fit)[1, 1])

naive <- coef(glm(hazard_window ~ exposed, data = d, family = binomial))[["exposed"]]
informative <- sum(tapply(d$exposed, d$patient, function(x) length(unique(x)) > 1))

z <- qnorm(0.975)
cat(sprintf("odds ratio %.2f (95%% CI %.2f to %.2f); %d of %d patients informative\n",
            exp(b), exp(b - z * se), exp(b + z * se), informative, length(unique(d$patient))))

cat("\n--- HARNESS ---\n")
cat(sprintf("log_or=%.10f\nlog_or_se=%.10f\nlog_or_lcl=%.10f\nlog_or_ucl=%.10f\n", b, se, b - z * se, b + z * se))
cat(sprintf("unconditional_log_or=%.10f\ninformative_patients=%d\nn_patients=%d\n", naive, informative,
            length(unique(d$patient))))
