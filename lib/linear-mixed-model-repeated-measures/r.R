# Linear mixed model for repeated measurements: the treatment-by-time interaction
#
# Each patient is measured several times, so their measurements are correlated. A random intercept per
# patient accounts for that; the treated-by-visit interaction is the difference in the rate of change
# between arms. A random intercept alone assumes every patient changes at the same rate apart from
# treatment; when individual trajectories differ, add a random slope and check the fit.

suppressPackageStartupMessages(library(lme4))

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d))

# PINNED: REML = TRUE, lmer's default and statsmodels' default, named because ML and REML give
# different variance components and therefore different standard errors. visit is NUMERIC (a linear
# trend); as a factor it would estimate a separate difference at every visit.
fit <- lmer(y ~ treated * visit + (1 | id), data = d, REML = TRUE)

co <- summary(fit)$coefficients
b <- co["treated:visit", "Estimate"]; se <- co["treated:visit", "Std. Error"]
vc <- as.data.frame(VarCorr(fit))
z <- qnorm(0.975)
cat(sprintf("difference in slope per visit %.3f, Wald 95%% %.3f to %.3f\n", b, b - z * se, b + z * se))
print(vc)

cat("\n--- HARNESS ---\n")
cat(sprintf("interaction=%.10f\ninteraction_se=%.10f\n", b, se))
cat(sprintf("interaction_lcl=%.10f\ninteraction_ucl=%.10f\n", b - z * se, b + z * se))
cat(sprintf("time_slope=%.10f\n", co["visit", "Estimate"]))
cat(sprintf("sd_patient=%.10f\nsd_residual=%.10f\n", vc$sdcor[vc$grp == "id"], vc$sdcor[vc$grp == "Residual"]))
cat(sprintf("n=%d\npatients=%d\n", nrow(d), length(unique(d$id))))
