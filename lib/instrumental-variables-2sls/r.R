# Instrumental variables: two-stage least squares (Mendelian-randomization style).
#
# 1. THE THREE ASSUMPTIONS ARE NOT TESTABLE FROM THESE DATA. The instrument must (i) move the exposure,
#    (ii) share no cause with the outcome, and (iii) affect the outcome only through the exposure. Only
#    (i) can be checked here, with the first-stage F statistic; (ii) and (iii) are arguments, and a
#    pleiotropic variant breaks (iii). State them in the protocol.
#
# 2. WEAK INSTRUMENTS BIAS 2SLS TOWARDS THE CONFOUNDED ESTIMATE. The first-stage F (ivreg's "Weak
#    instruments" diagnostic) should be well above 10; here it is about 300.
#
# 3. THE STANDARD ERROR FROM THE SECOND STAGE RUN BY HAND IS WRONG. Regressing Y on the fitted exposure
#    with lm() uses the wrong residuals (Y minus the fitted-exposure prediction, not Y minus the
#    exposure prediction). ivreg computes it correctly. Pinned: a heteroskedasticity-robust (HC1)
#    sandwich; the conventional SE is reported beside it.
#
# 4. THE ESTIMAND. With a constant effect, 2SLS estimates it; with effects that vary between people it
#    estimates a weighted average among those whose exposure the instrument moves.

library(ivreg)
library(sandwich)

d <- read.csv("fixture.csv")
fit <- ivreg(y ~ ldl | allele_score, data = d)
b <- coef(fit)[["ldl"]]
se_hc1 <- sqrt(vcovHC(fit, type = "HC1")["ldl", "ldl"])
se_conv <- sqrt(vcov(fit)["ldl", "ldl"])
diag_tab <- summary(fit, diagnostics = TRUE)$diagnostics
first_f <- diag_tab["Weak instruments", "statistic"]

ols <- coef(lm(y ~ ldl, data = d))[["ldl"]]

z <- qnorm(0.975)
cat(sprintf("2SLS effect %.3f (95%% CI %.3f to %.3f); OLS %.3f; first-stage F %.1f\n",
            b, b - z * se_hc1, b + z * se_hc1, ols, first_f))

cat("\n--- HARNESS ---\n")
cat(sprintf("iv_effect=%.10f\niv_se=%.10f\niv_lcl=%.10f\niv_ucl=%.10f\n", b, se_hc1, b - z * se_hc1, b + z * se_hc1))
cat(sprintf("iv_se_conventional=%.10f\n", se_conv))
cat(sprintf("first_stage_f=%.10f\nols_effect=%.10f\nn=%d\n", first_f, ols, nrow(d)))
