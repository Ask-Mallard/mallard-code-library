# Meta-analysis of diagnostic accuracy: the bivariate random-effects (Reitsma) model.
#
# 1. POOL SENSITIVITY AND SPECIFICITY JOINTLY. They are correlated across studies (a lower positivity
#    threshold raises sensitivity and lowers specificity), so pooling each separately, or pooling the
#    diagnostic odds ratio alone, loses information and can mislead. The bivariate model (Reitsma et al.
#    2005) treats the study-level logit sensitivity and logit specificity as bivariate normal, estimated
#    by REML (mada::reitsma).
#
# 2. REPORT the summary sensitivity and specificity with intervals, the between-study SDs and their
#    correlation, and a summary ROC curve or a prediction region; the between-study variation is usually
#    large and belongs in the conclusion.
#
# 3. Pinned: mada parameterizes (logit sensitivity, logit false positive rate); specificity is
#    1 - FPR. No zero cells here, so no continuity correction is applied.
#
# R ONLY: Python's standard stack has no bivariate random-effects meta-analysis; the control file refits
# the same model with metafor::rma.mv as an independent R implementation.

suppressPackageStartupMessages(library(mada))

d <- read.csv("fixture.csv")
fit <- reitsma(d[, c("TP", "FN", "FP", "TN")])
b <- coef(fit)
V <- vcov(fit)
psi <- fit$Psi
ls <- b[1, "tsens"]
lf <- b[1, "tfpr"]
se_ls <- sqrt(V[1, 1]); se_lf <- sqrt(V[2, 2])
z <- qnorm(0.975)

cat(sprintf("sensitivity %.3f (%.3f to %.3f); specificity %.3f (%.3f to %.3f); between-study SDs %.2f, %.2f; correlation %.2f\n",
            plogis(ls), plogis(ls - z * se_ls), plogis(ls + z * se_ls), 1 - plogis(lf), 1 - plogis(lf + z * se_lf),
            1 - plogis(lf - z * se_lf), sqrt(psi[1, 1]), sqrt(psi[2, 2]), -psi[1, 2] / sqrt(psi[1, 1] * psi[2, 2])))

cat("\n--- HARNESS ---\n")
cat(sprintf("logit_sensitivity=%.10f\nlogit_sensitivity_se=%.10f\n", ls, se_ls))
cat(sprintf("logit_specificity=%.10f\nlogit_specificity_se=%.10f\n", -lf, se_lf))
cat(sprintf("tau_sensitivity=%.10f\ntau_specificity=%.10f\n", sqrt(psi[1, 1]), sqrt(psi[2, 2])))
cat(sprintf("correlation_sens_spec=%.10f\n", -psi[1, 2] / sqrt(psi[1, 1] * psi[2, 2])))
cat(sprintf("studies=%d\n", nrow(d)))
