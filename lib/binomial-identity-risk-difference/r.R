# Adjusted risk difference from a binomial model with an identity link
#
# The treatment coefficient of an identity-link binomial model IS a risk difference, on the
# probability scale, adjusted for the covariates. The default (logit) link would give a log odds
# ratio instead, a different quantity.
#
# Identity-link models can fail to converge, or step outside [0, 1], when fitted risks approach 0 or 1.
# When that happens, estimate the risk difference by marginal standardization (g-computation) from
# a logistic model, or use a linear probability model with robust errors; do not drop the covariates
# to make it converge.

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d), all(d$outcome %in% c(0, 1)))

# PINNED: family = binomial(link = "identity") and start values. Without start, glm begins from the
# logit-scale default and often stops with "no valid set of coefficients has been found".
fit <- glm(outcome ~ treated + I(age - 60), family = binomial(link = "identity"), data = d,
           start = c(mean(d$outcome), 0, 0), control = glm.control(epsilon = 1e-12, maxit = 100))
stopifnot(fit$converged, all(fitted(fit) > 0 & fitted(fit) < 1))

est <- coef(summary(fit))
z <- qnorm(0.975)
rd <- est["treated", "Estimate"]; se <- est["treated", "Std. Error"]
print(est)

cat("\n--- HARNESS ---\n")
cat(sprintf("risk_difference=%.10f\nrisk_difference_se=%.10f\n", rd, se))
cat(sprintf("risk_difference_lcl=%.10f\nrisk_difference_ucl=%.10f\n", rd - z * se, rd + z * se))
cat(sprintf("age_slope=%.10f\n", est["I(age - 60)", "Estimate"]))
cat(sprintf("n=%d\n", nrow(d)))
