# Mann-Whitney (Wilcoxon rank-sum) test with the Hodges-Lehmann shift and its exact 95% interval
#
# The Hodges-Lehmann estimate is the median of all treated-minus-control pairwise differences. It is
# NOT the difference between the two medians, and the test does not compare medians unless the two
# groups have the same shape. Report the shift and its interval, not the p-value alone.

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d), all(d$group %in% c(0, 1)), !any(duplicated(d$y)))

treated <- d$y[d$group == 1]
control <- d$y[d$group == 0]

# PINNED: exact = TRUE and conf.int = TRUE. conf.int defaults to FALSE, so the estimate and interval
# are silently absent; exact defaults to TRUE only below 50 per group WITHOUT ties, and otherwise
# switches to a normal approximation with a continuity correction. Treated FIRST: the shift's sign
# follows argument order.
fit <- wilcox.test(treated, control, exact = TRUE, conf.int = TRUE, conf.level = 0.95)

cat(sprintf("W = %.0f, exact two-sided p = %.4g\n", fit$statistic, fit$p.value))
cat(sprintf("Hodges-Lehmann shift %.3f, exact 95%% %.3f to %.3f\n",
            fit$estimate, fit$conf.int[1], fit$conf.int[2]))

cat("\n--- HARNESS ---\n")
cat(sprintf("hl_shift=%.10f\n", fit$estimate))
cat(sprintf("hl_shift_lcl=%.10f\n", fit$conf.int[1]))
cat(sprintf("hl_shift_ucl=%.10f\n", fit$conf.int[2]))
cat(sprintf("u_statistic=%.10f\n", fit$statistic))
cat(sprintf("p_value=%.10e\n", fit$p.value))
cat(sprintf("n=%d\n", nrow(d)))
