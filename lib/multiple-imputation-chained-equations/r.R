# Multiple imputation by chained equations, pooled by Rubin's rules.
#
# 1. NAME THE MISSINGNESS ASSUMPTION FIRST. MI is valid under missing at random: missingness may depend
#    on OBSERVED variables. Here severity is missing more often when the outcome is high, which is MAR
#    given the outcome. Under missing not at random MI is not enough on its own; pair it with a
#    delta-adjusted sensitivity analysis (delta-adjusted-tipping-point).
#
# 2. THE OUTCOME GOES IN THE IMPUTATION MODEL. Leaving it out imputes the covariate as if unrelated to
#    the outcome and biases the coefficient towards zero. Complete-case analysis is biased here as well,
#    because the incomplete rows are not a random subset.
#
# 3. EVERY SETTING NAMED: m = 20 imputations (not mice's 5), maxit = 10 iterations, method "norm"
#    (Bayesian linear regression, which draws the parameters so the between-imputation variance is
#    right), and a seed. The analysis model is fitted on each completed dataset and pooled with Rubin's
#    rules; mice::pool gives Barnard-Rubin degrees of freedom.
#
# 4. REPORT the fraction of missing information beside the estimate.

library(mice)

d <- read.csv("fixture.csv")
imp <- mice(d[, c("treated", "severity", "y")], m = 20, maxit = 10, method = c("", "norm", ""),
            seed = 20261115, printFlag = FALSE)
fits <- with(imp, lm(y ~ treated + severity))
pooled <- pool(fits)
row <- pooled$pooled[pooled$pooled$term == "treated", ]
s <- summary(pooled, conf.int = TRUE)
srow <- s[s$term == "treated", ]

cc <- coef(lm(y ~ treated + severity, data = d))[["treated"]]

cat(sprintf("effect %.3f (95%% CI %.3f to %.3f), df %.1f, fraction of missing information %.2f; complete-case %.3f\n",
            row$estimate, srow$`2.5 %`, srow$`97.5 %`, row$df, row$fmi, cc))

cat("\n--- HARNESS ---\n")
cat(sprintf("effect=%.10f\neffect_se=%.10f\n", row$estimate, sqrt(row$t)))
cat(sprintf("effect_lcl=%.10f\neffect_ucl=%.10f\n", srow$`2.5 %`, srow$`97.5 %`))
cat(sprintf("within_variance=%.10f\nbetween_variance=%.10f\ndf=%.10f\nfmi=%.10f\n", row$ubar, row$b, row$df, row$fmi))
cat(sprintf("imputations=%d\ncomplete_case_effect=%.10f\n", row$m, cc))
cat(sprintf("missing=%d\nn=%d\n", sum(is.na(d$severity)), nrow(d)))
