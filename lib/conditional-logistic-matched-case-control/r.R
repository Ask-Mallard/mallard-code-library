# Conditional logistic regression for a 1:m matched case-control study
#
# Estimates the exposure odds ratio conditional on the matched sets. The matching factor is
# ELIMINATED from the likelihood rather than adjusted for, so it must not appear as a covariate:
# adding it would be conditioning on the thing already conditioned on.
#
# When the matched sets are risk sets (controls drawn from those still at risk at each case's event
# time, as in a nested case-control study), this conditional odds ratio estimates the cohort's hazard
# (rate) ratio, with no rare-disease assumption.

library(survival)

d <- read.csv("fixture.csv")

# PINNED DEFAULT: the reference level. R orders factor levels alphabetically, so a character
# exposure would silently pick its own reference and flip the odds ratio. Exposure is kept numeric
# 0/1 here, which is unambiguous; when your exposure is a factor, relevel() it explicitly.
d$exposed <- as.integer(d$exposed)

# PINNED DEFAULT: ties. clogit passes through to coxph, whose default is ties = "efron". For
# conditional logistic the EXACT method is the correct one and the two differ when a matched set
# has more than one case. Named here rather than inherited.
fit <- clogit(case ~ exposed + covariate + strata(set_id), data = d, method = "exact")

s <- summary(fit)
cat("estimates\n")
print(s$coefficients)

# Emit the machine-readable block the cross-language harness compares. Field names and order are
# fixed across all four languages so the comparison is positional-free.
cat("\n--- HARNESS ---\n")
cat(sprintf("exposure_log_or=%.10f\n", unname(coef(fit)["exposed"])))
cat(sprintf("exposure_se=%.10f\n", sqrt(diag(vcov(fit)))["exposed"]))
cat(sprintf("covariate_beta=%.10f\n", unname(coef(fit)["covariate"])))
cat(sprintf("n_sets=%d\n", length(unique(d$set_id))))
