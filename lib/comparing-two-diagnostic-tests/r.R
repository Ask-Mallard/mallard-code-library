# Comparing two diagnostic tests on the same patients: paired AUCs and paired sensitivities.
#
# 1. PAIRED DATA NEED A PAIRED COMPARISON. Both tests are measured on everyone, so their errors are
#    correlated. DeLong's test for two correlated ROC curves uses the covariance of the two AUCs; treating
#    them as independent overstates the standard error of the difference when the markers are positively
#    correlated (reported beside it).
#
# 2. AT A THRESHOLD, compare sensitivities among the patients WITH disease, and specificities among those
#    without, by McNemar's test on the discordant pairs. A two-proportion test ignores the pairing. The
#    thresholds must be prespecified; choosing each test's best threshold on these data favours both.
#    Pinned: McNemar without continuity correction; the paired difference in sensitivity with a Wald
#    interval for paired proportions.
#
# 3. The direction is pinned for both markers: higher values mean disease.

suppressPackageStartupMessages(library(pROC))

d <- read.csv("fixture.csv")
ra <- roc(d$disease, d$marker_a, levels = c(0, 1), direction = "<", quiet = TRUE)
rb <- roc(d$disease, d$marker_b, levels = c(0, 1), direction = "<", quiet = TRUE)
test <- roc.test(ra, rb, method = "delong", paired = TRUE)
diff <- as.numeric(ra$auc) - as.numeric(rb$auc)
se_paired <- sqrt(var(ra, method = "delong") + var(rb, method = "delong") - 2 * cov(ra, rb, method = "delong"))
se_unpaired <- sqrt(var(ra, method = "delong") + var(rb, method = "delong"))

# Sensitivities at the prespecified thresholds, among patients with disease.
cases <- d[d$disease == 1, ]
pa <- as.integer(cases$marker_a > 0.6)
pb <- as.integer(cases$marker_b > 0.6)
tab <- table(factor(pa, 0:1), factor(pb, 0:1))
mc <- mcnemar.test(tab, correct = FALSE)
n <- nrow(cases)
b01 <- tab["1", "0"]; b10 <- tab["0", "1"]
sens_diff <- (b01 - b10) / n
sens_diff_se <- sqrt((b01 + b10) - (b01 - b10)^2 / n) / n

z <- qnorm(0.975)
cat(sprintf("AUC A %.3f, B %.3f, difference %.3f (paired DeLong p = %.4f); sensitivity difference %.3f (McNemar p = %.4f)\n",
            as.numeric(ra$auc), as.numeric(rb$auc), diff, test$p.value, sens_diff, mc$p.value))

cat("\n--- HARNESS ---\n")
cat(sprintf("auc_a=%.10f\nauc_b=%.10f\nauc_difference=%.10f\n", as.numeric(ra$auc), as.numeric(rb$auc), diff))
cat(sprintf("auc_difference_se=%.10f\nauc_difference_se_unpaired=%.10f\nauc_test_z=%.10f\n", se_paired, se_unpaired, test$statistic))
cat(sprintf("sensitivity_a=%.10f\nsensitivity_b=%.10f\nsensitivity_difference=%.10f\nsensitivity_difference_se=%.10f\n",
            mean(pa), mean(pb), sens_diff, sens_diff_se))
cat(sprintf("mcnemar_statistic=%.10f\ndiscordant_a_only=%d\ndiscordant_b_only=%d\nn=%d\n", mc$statistic, b01, b10, nrow(d)))
