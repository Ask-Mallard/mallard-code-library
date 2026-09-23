# McNemar's test for a paired binary outcome, with the change in proportion and its 95% interval
#
# The same people tested twice. Only the discordant pairs (positive then negative, negative then
# positive) carry information about change. The estimate is the change in the proportion positive,
# (c - b) / n, with the paired Wald interval; the tests are McNemar's chi-square and the exact
# binomial test on the discordant pairs.

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d), all(d$before %in% c(0, 1)), all(d$after %in% c(0, 1)))

n <- nrow(d)
b <- sum(d$before == 1 & d$after == 0)   # positive -> negative
c <- sum(d$before == 0 & d$after == 1)   # negative -> positive
tab <- matrix(c(sum(d$before == 1 & d$after == 1), b, c, sum(d$before == 0 & d$after == 0)),
              nrow = 2, byrow = TRUE)

# PINNED: correct = FALSE. mcnemar.test applies a continuity correction by default; the exact
# binomial p-value below is what a small discordant count calls for instead.
chi <- mcnemar.test(tab, correct = FALSE)
exact <- binom.test(b, b + c, p = 0.5)

change <- (c - b) / n
se <- sqrt(b + c - (c - b)^2 / n) / n
z <- qnorm(0.975)

cat(sprintf("discordant: %d positive->negative, %d negative->positive\n", b, c))
cat(sprintf("change in proportion positive %.4f, 95%% %.4f to %.4f\n", change, change - z * se, change + z * se))
cat(sprintf("McNemar chi-square %.3f; exact p %.4g\n", chi$statistic, exact$p.value))

cat("\n--- HARNESS ---\n")
cat(sprintf("change=%.10f\nchange_lcl=%.10f\nchange_ucl=%.10f\n", change, change - z * se, change + z * se))
cat(sprintf("discordant_pos_neg=%d\ndiscordant_neg_pos=%d\n", b, c))
cat(sprintf("mcnemar_chi2=%.10f\n", chi$statistic))
cat(sprintf("exact_p_value=%.10e\n", exact$p.value))
cat(sprintf("n=%d\n", n))
