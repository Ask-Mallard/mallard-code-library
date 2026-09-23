# Self-controlled case series (SCCS): conditional Poisson regression within each case.
#
# 1. ONLY CASES, EACH THEIR OWN CONTROL. The event rate in the risk window after exposure is compared
#    with the rate at other times in the same person's observation period. Everything fixed over that
#    period (sex, genes, underlying frailty) is removed, measured or not.
#
# 2. THE MODEL CONDITIONS ON EACH PERSON'S NUMBER OF EVENTS. That conditional (multinomial) likelihood
#    has the same estimate and standard error as a Poisson regression with one fixed effect per person
#    and log(interval length) as the offset, which is what this file fits (the control file fits the
#    multinomial form directly and checks it). The SCCS package's standardsccs() does the same.
#
# 3. AGE (OR CALENDAR TIME) MUST BE IN THE MODEL. Exposure is usually given at particular ages, and event
#    rates change with age. Here vaccination falls in the high-risk age band, and leaving age out
#    credits the age effect to the vaccine.
#
# 4. ASSUMPTIONS: events do not change the probability of later exposure (or the observation period),
#    and recurrent events are independent, or the first event only is used for a non-recurrent one.

d <- read.csv("fixture.csv")
fit <- glm(events ~ risk + factor(age_band) + factor(child) + offset(log(days)), data = d, family = poisson)
b <- coef(fit)[["risk"]]
se <- sqrt(vcov(fit)["risk", "risk"])

no_age <- coef(glm(events ~ risk + factor(child) + offset(log(days)), data = d, family = poisson))[["risk"]]

z <- qnorm(0.975)
cat(sprintf("incidence rate ratio %.2f (95%% CI %.2f to %.2f)\n", exp(b), exp(b - z * se), exp(b + z * se)))

cat("\n--- HARNESS ---\n")
cat(sprintf("log_irr=%.10f\nlog_irr_se=%.10f\nlog_irr_lcl=%.10f\nlog_irr_ucl=%.10f\n", b, se, b - z * se, b + z * se))
cat(sprintf("log_irr_age_ignored=%.10f\n", no_age))
cat(sprintf("cases=%d\nevents=%d\n", length(unique(d$child)), sum(d$events)))
