# Predictive values and likelihood ratios from a diagnostic accuracy study.
#
# 1. PREDICTIVE VALUES DEPEND ON PREVALENCE. PPV and NPV estimated from a cohort apply to a population with
#    that cohort's prevalence. For another setting, carry the test's sensitivity and specificity (or its
#    likelihood ratios) across and recompute: PPV = se p / (se p + (1 - sp)(1 - p)). A case-control
#    sample, whose "prevalence" is set by design, gives no valid predictive values at all.
#
# 2. INTERVALS MATCH THE ESTIMATOR. PPV and NPV are single proportions among the test-positives and
#    test-negatives: Wilson score intervals (prop.test without continuity correction). Likelihood ratios
#    are ratios of two independent proportions: a log-scale interval, SE of log LR+ =
#    sqrt((1 - se)/(se n_d) + sp/((1 - sp) n_nd)) (Simel, Samsa and Matchar 1991).
#
# 3. The PPV at a stated other prevalence is reported to make point 1 concrete.

d <- read.csv("fixture.csv")
tp <- sum(d$disease == 1 & d$test_positive == 1)
fn <- sum(d$disease == 1 & d$test_positive == 0)
fp <- sum(d$disease == 0 & d$test_positive == 1)
tn <- sum(d$disease == 0 & d$test_positive == 0)
nd <- tp + fn
nnd <- fp + tn

wilson <- function(x, n) as.numeric(prop.test(x, n, correct = FALSE)$conf.int)
ppv <- tp / (tp + fp); ppv_ci <- wilson(tp, tp + fp)
npv <- tn / (tn + fn); npv_ci <- wilson(tn, tn + fn)
se <- tp / nd
sp <- tn / nnd
lr_pos <- se / (1 - sp)
lr_neg <- (1 - se) / sp
se_log_pos <- sqrt((1 - se) / (se * nd) + sp / ((1 - sp) * nnd))
se_log_neg <- sqrt(se / ((1 - se) * nd) + (1 - sp) / (sp * nnd))
p_other <- 0.05
ppv_other <- se * p_other / (se * p_other + (1 - sp) * (1 - p_other))

z <- qnorm(0.975)
cat(sprintf("PPV %.3f (%.3f to %.3f), NPV %.3f (%.3f to %.3f); LR+ %.2f (%.2f to %.2f), LR- %.3f (%.3f to %.3f); PPV at 5%% prevalence %.3f\n",
            ppv, ppv_ci[1], ppv_ci[2], npv, npv_ci[1], npv_ci[2], lr_pos, exp(log(lr_pos) - z * se_log_pos),
            exp(log(lr_pos) + z * se_log_pos), lr_neg, exp(log(lr_neg) - z * se_log_neg), exp(log(lr_neg) + z * se_log_neg), ppv_other))

cat("\n--- HARNESS ---\n")
cat(sprintf("ppv=%.10f\nppv_lcl=%.10f\nppv_ucl=%.10f\n", ppv, ppv_ci[1], ppv_ci[2]))
cat(sprintf("npv=%.10f\nnpv_lcl=%.10f\nnpv_ucl=%.10f\n", npv, npv_ci[1], npv_ci[2]))
cat(sprintf("log_lr_positive=%.10f\nlog_lr_positive_se=%.10f\n", log(lr_pos), se_log_pos))
cat(sprintf("log_lr_negative=%.10f\nlog_lr_negative_se=%.10f\n", log(lr_neg), se_log_neg))
cat(sprintf("ppv_at_prevalence_0.05=%.10f\n", ppv_other))
cat(sprintf("n=%d\n", nrow(d)))
