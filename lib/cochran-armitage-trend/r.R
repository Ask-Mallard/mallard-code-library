# Cochran-Armitage test for trend in a proportion across ordered groups, with the linear slope
#
# Asks whether the outcome proportion rises (or falls) steadily with an ordered exposure, on the
# scores given. It uses 1 degree of freedom where a chi-square test of the whole table uses k - 1,
# so it is more powerful when the trend is real, and blind to a non-monotone pattern. The slope is
# the change in proportion per score step, from a weighted linear fit.

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d), all(d$outcome %in% c(0, 1)))

score <- sort(unique(d$dose))
events <- as.vector(tapply(d$outcome, d$dose, sum)[as.character(score)])
totals <- as.vector(table(d$dose)[as.character(score)])

# PINNED: score = the dose values themselves. prop.trend.test's default score is seq_along(x),
# 1, 2, ..., k, which matches equally spaced doses only by coincidence of spacing.
trend <- prop.trend.test(events, totals, score = score)
slope <- coef(lm(I(events / totals) ~ score, weights = totals))[["score"]]

print(data.frame(score, events, totals, proportion = events / totals))
cat(sprintf("trend chi-square %.4f on 1 df, p = %.4g; slope %.4f per step\n",
            trend$statistic, trend$p.value, slope))

cat("\n--- HARNESS ---\n")
cat(sprintf("trend_chi2=%.10f\n", trend$statistic))
cat(sprintf("p_value=%.10e\n", trend$p.value))
cat(sprintf("slope=%.10f\n", slope))
cat(sprintf("n=%d\n", nrow(d)))
