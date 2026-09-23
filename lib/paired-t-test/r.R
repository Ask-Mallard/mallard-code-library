# Paired t-test: the mean within-person change, after minus before, with a 95% interval
#
# Each person is their own control, so the analysis is a one-sample t on the differences. A
# two-sample (Welch) test on the same two columns ignores the pairing, treats the between-person
# spread as noise, and gives a far wider interval for the same data.
#
# Assumes the DIFFERENCES are roughly normal (not the raw values). With a skewed or outlying set of
# differences, see wilcoxon-signed-rank-hodges-lehmann.

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d))

# PINNED: paired = TRUE, with after FIRST so the difference is after minus before. The sign of the
# reported change follows argument order, and a reversed order silently flips "fell by 5" into
# "rose by 5". R 4.4 also removed paired = TRUE from the formula interface, so the two-vector form
# is the one that still works.
fit <- t.test(d$after, d$before, paired = TRUE, conf.level = 0.95)

cat(sprintf("%d pairs; mean change (after - before) %.2f, 95%% %.2f to %.2f; t = %.3f on %d df\n",
            nrow(d), fit$estimate, fit$conf.int[1], fit$conf.int[2], fit$statistic, as.integer(fit$parameter)))

cat("\n--- HARNESS ---\n")
cat(sprintf("mean_change=%.10f\n", fit$estimate))
cat(sprintf("mean_change_lcl=%.10f\n", fit$conf.int[1]))
cat(sprintf("mean_change_ucl=%.10f\n", fit$conf.int[2]))
cat(sprintf("sd_change=%.10f\n", sd(d$after - d$before)))
cat(sprintf("t_statistic=%.10f\n", fit$statistic))
cat(sprintf("df=%d\n", as.integer(fit$parameter)))
cat(sprintf("n=%d\n", nrow(d)))
