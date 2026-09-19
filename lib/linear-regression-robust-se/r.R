# Linear regression with heteroskedasticity-consistent standard errors
#
# The coefficients are ordinary least squares and are unaffected by everything below. What the
# robust variance changes is the STANDARD ERROR, and therefore the confidence interval and the
# p value. Unequal spread between arms does not bias the estimate; it makes the classical interval
# the wrong width.

library(sandwich)
library(lmtest)

d <- read.csv("fixture.csv")

fit <- lm(y ~ exposed + covariate, data = d)

# PINNED DEFAULT: type = "HC3".
#
# Every package has a button labelled "robust" and each one presses a DIFFERENT estimator:
#   R    sandwich::vcovHC   defaults to HC3
#   Python statsmodels      has no default at all -- .fit() alone is the CLASSICAL variance
#
# So the two files disagree by default -- R gives HC3, a bare statsmodels .fit() the classical
# variance -- and nothing in either output says which one it used. Every file in this entry names HC3.
V_hc3 <- vcovHC(fit, type = "HC3")

# The ladder, printed so a reader can see the size of what the choice costs rather than take it
# on trust. On this fixture the classical SE for the exposure is about 7% smaller than HC3, while
# HC0, HC1 and HC3 sit within 2e-4 of one another.
V_hc0 <- vcovHC(fit, type = "HC0")
V_hc1 <- vcovHC(fit, type = "HC1")

cat("estimates with HC3 standard errors\n")
print(coeftest(fit, vcov. = V_hc3))
cat("\nclassical standard errors, for comparison ONLY -- they assume constant variance\n")
print(summary(fit)$coefficients)

cat("\n--- HARNESS ---\n")
cat(sprintf("intercept=%.10f\n", unname(coef(fit)["(Intercept)"])))
cat(sprintf("exposure_beta=%.10f\n", unname(coef(fit)["exposed"])))
cat(sprintf("covariate_beta=%.10f\n", unname(coef(fit)["covariate"])))
cat(sprintf("exposure_se=%.10f\n", sqrt(V_hc3["exposed", "exposed"])))
cat(sprintf("covariate_se=%.10f\n", sqrt(V_hc3["covariate", "covariate"])))
cat(sprintf("exposure_se_hc0=%.10f\n", sqrt(V_hc0["exposed", "exposed"])))
cat(sprintf("exposure_se_hc1=%.10f\n", sqrt(V_hc1["exposed", "exposed"])))
cat(sprintf("exposure_se_classical=%.10f\n", summary(fit)$coefficients["exposed", "Std. Error"]))
cat(sprintf("n=%d\n", nrow(d)))
