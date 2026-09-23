# Stepped-wedge trial: the Hussey-Hughes linear mixed model
#
# Period as FIXED effects and a random intercept per ward. The period effects are what separate the
# intervention from the secular trend: in a stepped wedge the intervention periods are the later ones,
# so without them any improvement over time is credited to the intervention.
#
# With few wards (here 12) the Wald interval below is too narrow; a Kenward-Roger or Satterthwaite
# correction (lmerTest), or randomization inference over the wedge's possible orders, is the honest
# small-sample report. This entry's Wald interval is the cross-language baseline.

suppressPackageStartupMessages(library(lme4))

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d))

# PINNED: factor(period) (one effect per period, not a straight line) and REML = TRUE, named.
fit <- lmer(y ~ intervention + factor(period) + (1 | ward), data = d, REML = TRUE)

co <- summary(fit)$coefficients
b <- co["intervention", "Estimate"]; se <- co["intervention", "Std. Error"]
vc <- as.data.frame(VarCorr(fit))
z <- qnorm(0.975)
cat(sprintf("intervention effect %.3f, Wald 95%% %.3f to %.3f\n", b, b - z * se, b + z * se))

cat("\n--- HARNESS ---\n")
cat(sprintf("effect=%.10f\neffect_se=%.10f\neffect_lcl=%.10f\neffect_ucl=%.10f\n", b, se, b - z * se, b + z * se))
cat(sprintf("sd_ward=%.10f\nsd_residual=%.10f\n", vc$sdcor[vc$grp == "ward"], vc$sdcor[vc$grp == "Residual"]))
cat(sprintf("n=%d\nwards=%d\n", nrow(d), length(unique(d$ward))))
