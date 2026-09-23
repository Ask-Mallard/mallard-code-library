# ANCOVA: the treatment effect on a follow-up measurement, adjusted for its baseline value
#
# In a randomized trial with a baseline measurement, regress the follow-up value on treatment and
# baseline. It is unbiased like the alternatives (follow-up only, or change score) and more precise
# than either, because it estimates how much of the baseline carries forward instead of assuming 0
# or 1.
#
# Adjust only for what was measured BEFORE randomization. A variable measured after treatment starts
# can be affected by it, and adjusting for it biases the effect.
#
# Model-based (classical) standard errors, appropriate with equal allocation and similar spread in
# both arms. With clearly unequal variances or unequal allocation, use heteroskedasticity-robust
# errors (see linear-regression-robust-se).

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d), all(d$treated %in% c(0, 1)))

# PINNED: the outcome is the FOLLOW-UP value with baseline as a covariate. Using the change score as
# the outcome AND adjusting for baseline gives the same treatment estimate, but a change score
# WITHOUT baseline adjustment forces the carry-over to 1 and loses precision.
fit <- lm(followup ~ treated + baseline, data = d)
ci <- confint(fit, level = 0.95)

print(summary(fit)$coefficients)

cat("\n--- HARNESS ---\n")
cat(sprintf("effect=%.10f\n", coef(fit)[["treated"]]))
cat(sprintf("effect_se=%.10f\n", summary(fit)$coefficients["treated", "Std. Error"]))
cat(sprintf("effect_lcl=%.10f\neffect_ucl=%.10f\n", ci["treated", 1], ci["treated", 2]))
cat(sprintf("baseline_coef=%.10f\n", coef(fit)[["baseline"]]))
cat(sprintf("n=%d\n", nrow(d)))
