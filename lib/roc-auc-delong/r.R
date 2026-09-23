# ROC analysis: the area under the curve with a DeLong interval, and a threshold.
#
# 1. THE AUC is the probability that a randomly chosen patient with disease has a higher marker value
#    than one without. It summarizes discrimination over all thresholds; it says nothing about
#    calibration or about how the test performs at the threshold a clinic would use.
#
# 2. NAME THE DIRECTION. pROC's default (direction = "auto") picks whichever direction gives AUC > 0.5,
#    which silently turns a marker that points the wrong way into a good one. Pinned: higher values mean
#    disease (direction = "<", controls first).
#
# 3. THE INTERVAL is DeLong's nonparametric one (ci.auc method = "delong"), on the AUC scale.
#
# 4. THE THRESHOLD maximizing Youden's index is reported with its sensitivity and specificity. It is
#    chosen on these data, so its performance here is optimistic; a threshold for use should be
#    prespecified or validated elsewhere. pROC's candidate thresholds are midpoints between consecutive
#    observed values.

suppressPackageStartupMessages(library(pROC))

d <- read.csv("fixture.csv")
r <- roc(d$disease, d$marker, levels = c(0, 1), direction = "<", quiet = TRUE)
ci <- ci.auc(r, method = "delong")
se <- sqrt(var(r, method = "delong"))
best <- coords(r, "best", best.method = "youden", ret = c("threshold", "sensitivity", "specificity"))

cat(sprintf("AUC %.3f (DeLong 95%% CI %.3f to %.3f); Youden threshold %.3f (sensitivity %.3f, specificity %.3f)\n",
            as.numeric(r$auc), ci[1], ci[3], best$threshold, best$sensitivity, best$specificity))

cat("\n--- HARNESS ---\n")
cat(sprintf("auc=%.10f\nauc_se=%.10f\nauc_lcl=%.10f\nauc_ucl=%.10f\n", as.numeric(r$auc), se, ci[1], ci[3]))
cat(sprintf("youden_threshold=%.10f\nsensitivity=%.10f\nspecificity=%.10f\n", best$threshold, best$sensitivity, best$specificity))
cat(sprintf("cases=%d\ncontrols=%d\n", sum(d$disease == 1), sum(d$disease == 0)))
