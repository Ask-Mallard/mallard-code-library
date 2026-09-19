# Two-by-two table: chi-square, Fisher's exact test and the odds ratio
#
# Base R only. No package is needed for any of this, which is worth noticing: what differs between
# languages here is not capability but which test each one runs when you do not say.

d <- read.csv("fixture.csv")

# Levels are stated explicitly so the table is [exposed, unexposed] x [outcome, no outcome] rather
# than whatever order the factor happens to take. read.csv gives integers, which sort 0 then 1, so
# the table would otherwise be upside down and every quantity below would refer to the wrong cell
# while looking entirely reasonable.
tbl <- table(factor(d$exposed, levels = c(1, 0)),
             factor(d$outcome, levels = c(1, 0)))
a <- tbl[1, 1]; b <- tbl[1, 2]; c_ <- tbl[2, 1]; dd <- tbl[2, 2]

# PINNED DEFAULT: correct = FALSE, the uncorrected Pearson chi-square.
#
# R and scipy BOTH apply Yates' continuity correction to a 2x2 by default, so both files switch it
# off to report the uncorrected Pearson. On this fixture that is 49.90 uncorrected against 49.19
# corrected -- with every expected count in the hundreds, which is not the sparse-table edge case
# the correction is usually defended for.
#
# Uncorrected is pinned because the correction is a small-sample device for approximating an exact
# conditional test, and it is conservative to a fault where the expected counts are large. Where
# they are NOT large, the honest answer is Fisher's exact test, computed below, rather than a
# corrected approximation.
pearson <- chisq.test(tbl, correct = FALSE)
yates <- chisq.test(tbl, correct = TRUE)
fisher <- fisher.test(tbl)

cat("2x2 table (rows: exposed, unexposed; columns: outcome, no outcome)\n")
print(tbl)
cat("\nPearson chi-square, uncorrected\n"); print(pearson)
cat("\nwith Yates' continuity correction, FOR COMPARISON ONLY\n"); print(yates)
cat("\nFisher's exact test\n"); print(fisher)

# THE ODDS RATIO IS COMPUTED FROM THE TABLE, NOT TAKEN FROM fisher.test.
# fisher.test()$estimate is the CONDITIONAL MAXIMUM LIKELIHOOD odds ratio, which is a different
# estimator from the sample odds ratio ad/bc, and scipy's fisher_exact returns the sample one. Two
# files each reporting "the odds ratio from Fisher's test" would therefore disagree, and the
# disagreement would look like a bug in one of them rather than two different estimators.
log_or <- log((a * dd) / (b * c_))
log_or_se <- sqrt(1 / a + 1 / b + 1 / c_ + 1 / dd)
cat(sprintf("\nsample odds ratio %.6f (95%% CI %.6f to %.6f)\n", exp(log_or),
            exp(log_or - 1.96 * log_or_se), exp(log_or + 1.96 * log_or_se)))
cat(sprintf("conditional MLE odds ratio from fisher.test: %.6f -- A DIFFERENT ESTIMATOR\n",
            unname(fisher$estimate)))

cat("\n--- HARNESS ---\n")
cat(sprintf("chisq_uncorrected=%.10f\n", unname(pearson$statistic)))
cat(sprintf("chisq_yates=%.10f\n", unname(yates$statistic)))
cat(sprintf("log_or=%.10f\n", log_or))
cat(sprintf("log_or_se=%.10f\n", log_or_se))
cat(sprintf("neg_log10_p_fisher=%.10f\n", -log10(fisher$p.value)))
cat(sprintf("risk_exposed=%.10f\n", a / (a + b)))
cat(sprintf("risk_unexposed=%.10f\n", c_ / (c_ + dd)))
cat(sprintf("risk_difference=%.10f\n", a / (a + b) - c_ / (c_ + dd)))
cat(sprintf("n=%d\n", nrow(d)))
