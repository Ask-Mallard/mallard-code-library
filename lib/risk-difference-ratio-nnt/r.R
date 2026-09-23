# Risk difference, risk ratio and number needed to treat from a two-arm trial, with 95% intervals
#
# Absolute and relative effects together: the risk difference (treated minus control) with Newcombe's
# hybrid score interval, the risk ratio with the log (Katz) interval, and the NNT with the interval
# obtained by inverting the risk difference's bounds. A relative effect alone hides the baseline risk
# that decides whether an effect matters to a patient.

library(binom)

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d), all(d$treated %in% c(0, 1)), all(d$outcome %in% c(0, 1)))

x1 <- sum(d$outcome[d$treated == 1]); n1 <- sum(d$treated == 1)
x0 <- sum(d$outcome[d$treated == 0]); n0 <- sum(d$treated == 0)
p1 <- x1 / n1; p0 <- x0 / n0
z <- qnorm(0.975)

# Newcombe (1998) method 10: combine each arm's Wilson interval. PINNED: methods = "wilson".
w1 <- binom.confint(x1, n1, methods = "wilson"); w0 <- binom.confint(x0, n0, methods = "wilson")
rd <- p1 - p0
rd_lcl <- rd - sqrt((p1 - w1$lower)^2 + (w0$upper - p0)^2)
rd_ucl <- rd + sqrt((w1$upper - p1)^2 + (p0 - w0$lower)^2)

# Risk ratio on the log scale (Katz).
log_rr <- log(p1 / p0)
se_log_rr <- sqrt(1 / x1 - 1 / n1 + 1 / x0 - 1 / n0)

# NNT = 1 / |RD|. Its interval inverts the RD bounds ONLY when they exclude zero; otherwise it runs
# through infinity and must be reported as two ranges, never as a single interval.
stopifnot(rd_ucl < 0 || rd_lcl > 0)
nnt <- 1 / abs(rd)
nnt_bounds <- sort(1 / abs(c(rd_lcl, rd_ucl)))

cat(sprintf("treated %d/%d (%.3f), control %d/%d (%.3f)\n", x1, n1, p1, x0, n0, p0))
cat(sprintf("risk difference %.4f, Newcombe 95%% %.4f to %.4f\n", rd, rd_lcl, rd_ucl))
cat(sprintf("risk ratio %.3f, 95%% %.3f to %.3f\n", exp(log_rr), exp(log_rr - z * se_log_rr), exp(log_rr + z * se_log_rr)))
cat(sprintf("NNT %.1f, 95%% %.1f to %.1f\n", nnt, nnt_bounds[1], nnt_bounds[2]))

cat("\n--- HARNESS ---\n")
cat(sprintf("risk_treated=%.10f\nrisk_control=%.10f\n", p1, p0))
cat(sprintf("risk_difference=%.10f\nrisk_difference_lcl=%.10f\nrisk_difference_ucl=%.10f\n", rd, rd_lcl, rd_ucl))
cat(sprintf("log_risk_ratio=%.10f\n", log_rr))
cat(sprintf("risk_ratio=%.10f\nrisk_ratio_lcl=%.10f\nrisk_ratio_ucl=%.10f\n",
            exp(log_rr), exp(log_rr - z * se_log_rr), exp(log_rr + z * se_log_rr)))
cat(sprintf("nnt=%.10f\nnnt_lcl=%.10f\nnnt_ucl=%.10f\n", nnt, nnt_bounds[1], nnt_bounds[2]))
cat(sprintf("n=%d\n", nrow(d)))
