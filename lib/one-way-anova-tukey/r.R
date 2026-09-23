# One-way ANOVA across three groups, then all pairwise differences with Tukey-adjusted 95% intervals
#
# The F test asks whether ANY group mean differs; it does not say which. The pairwise differences
# are the result a reader needs, and with three comparisons their intervals are adjusted (Tukey's
# honestly significant difference) so the family of three keeps 95% coverage together.
#
# Assumes roughly normal outcomes with similar SDs across groups. With clearly unequal SDs, prefer
# Welch's ANOVA (oneway.test) and Games-Howell comparisons; with skewed outcomes, see kruskal-wallis.

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d))

# PINNED: factor(arm). arm is coded 1, 2, 3, and without factor() aov fits a single numeric slope
# on 1 df, reported as an ANOVA table that looks entirely normal.
fit <- aov(y ~ factor(arm), data = d)
tab <- summary(fit)[[1]]
tk <- TukeyHSD(fit, conf.level = 0.95)[[1]]

print(tab)
print(tk)

cat("\n--- HARNESS ---\n")
cat(sprintf("f_statistic=%.10f\n", tab[1, "F value"]))
cat(sprintf("df_between=%d\ndf_within=%d\n", as.integer(tab[1, "Df"]), as.integer(tab[2, "Df"])))
cat(sprintf("pooled_sd=%.10f\n", sqrt(tab[2, "Mean Sq"])))
for (pair in c("2-1", "3-1", "3-2")) {
  key <- paste0("diff_", sub("-", "_", pair))
  cat(sprintf("%s=%.10f\n%s_lcl=%.10f\n%s_ucl=%.10f\n",
              key, tk[pair, "diff"], key, tk[pair, "lwr"], key, tk[pair, "upr"]))
}
cat(sprintf("n=%d\n", nrow(d)))
