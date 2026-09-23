# Small-sample cluster-robust inference: CR2 standard errors with Satterthwaite degrees of freedom
#
# A patient-level regression in a trial with few clusters. The ordinary cluster-robust sandwich (CR0 or
# CR1) is too small with few clusters, and a normal reference makes it worse. CR2, the Bell-McCaffrey
# bias-reduced sandwich, with Satterthwaite degrees of freedom, keeps the interval close to its stated
# coverage (Pustejovsky and Tipton 2018).

suppressPackageStartupMessages(library(clubSandwich))

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d))

fit <- lm(y ~ treated + I(age - 50), data = d)

# PINNED: vcov = "CR2" and test = "Satterthwaite". clubSandwich also offers CR0 and CR1 (the Stata-style
# small-sample factor), which are smaller here, and a naive t or z reference.
ct <- coef_test(fit, vcov = "CR2", cluster = d$clinic, test = "Satterthwaite")
ci <- conf_int(fit, vcov = "CR2", cluster = d$clinic, test = "Satterthwaite", level = 0.95)
row <- which(ct$Coef == "treated")

cat(sprintf("treatment effect %.3f, CR2 SE %.4f, Satterthwaite df %.2f, 95%% %.3f to %.3f\n",
            ct$beta[row], ct$SE[row], ct$df_Satt[row], ci$CI_L[row], ci$CI_U[row]))

cat("\n--- HARNESS ---\n")
cat(sprintf("effect=%.10f\neffect_cr2_se=%.10f\neffect_df=%.10f\n", ct$beta[row], ct$SE[row], ct$df_Satt[row]))
cat(sprintf("effect_lcl=%.10f\neffect_ucl=%.10f\n", ci$CI_L[row], ci$CI_U[row]))
cat(sprintf("clusters=%d\nn=%d\n", length(unique(d$clinic)), nrow(d)))
