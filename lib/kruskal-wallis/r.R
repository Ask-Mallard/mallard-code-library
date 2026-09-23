# Kruskal-Wallis test across three groups, reported with each group's median
#
# The rank-based counterpart of one-way ANOVA. The test says whether the groups differ in
# distribution; it gives no effect size, so the group medians are reported beside it and pairwise
# differences, if needed, come from mann-whitney-hodges-lehmann with a multiplicity adjustment.
# It reads as a comparison of medians only when the groups share one shape.

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d))

# PINNED: factor(arm). kruskal.test's formula interface groups by the factor's levels; a numeric
# arm column is accepted too, but pinning the factor keeps the grouping explicit and matches the
# ANOVA entry. The chi-square approximation (with the tie correction) is the default and is used.
fit <- kruskal.test(y ~ factor(arm), data = d)
med <- tapply(d$y, d$arm, median)

print(fit)
print(med)

cat("\n--- HARNESS ---\n")
cat(sprintf("h_statistic=%.10f\n", fit$statistic))
cat(sprintf("df=%d\n", as.integer(fit$parameter)))
cat(sprintf("p_value=%.10e\n", fit$p.value))
cat(sprintf("median_arm1=%.10f\nmedian_arm2=%.10f\nmedian_arm3=%.10f\n", med["1"], med["2"], med["3"]))
cat(sprintf("n=%d\n", nrow(d)))
