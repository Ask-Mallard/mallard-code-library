# Wilcoxon signed-rank test with the Hodges-Lehmann pseudo-median change and its exact 95% interval
#
# For paired data whose changes are symmetric but not normal (heavy tails, a few large changes).
# The pseudo-median is the median of all pairwise averages of the changes (Walsh averages); under
# symmetry it estimates the centre of the change distribution. If the changes are clearly skewed,
# the test no longer answers a question about a typical change, and neither does the estimate.

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d))

# PINNED: paired = TRUE, exact = TRUE, conf.int = TRUE, after FIRST. conf.int defaults to FALSE;
# exact falls back to a normal approximation with zeros, ties or 50 or more pairs; and the sign of the
# change follows argument order.
fit <- wilcox.test(d$after, d$before, paired = TRUE, exact = TRUE, conf.int = TRUE, conf.level = 0.95)

cat(sprintf("V = %.0f, exact two-sided p = %.4g\n", fit$statistic, fit$p.value))
cat(sprintf("pseudo-median change %.3f, exact 95%% %.3f to %.3f\n",
            fit$estimate, fit$conf.int[1], fit$conf.int[2]))

cat("\n--- HARNESS ---\n")
cat(sprintf("pseudomedian_change=%.10f\n", fit$estimate))
cat(sprintf("pseudomedian_change_lcl=%.10f\n", fit$conf.int[1]))
cat(sprintf("pseudomedian_change_ucl=%.10f\n", fit$conf.int[2]))
cat(sprintf("v_statistic=%.10f\n", fit$statistic))
cat(sprintf("p_value=%.10e\n", fit$p.value))
cat(sprintf("n=%d\n", nrow(d)))
