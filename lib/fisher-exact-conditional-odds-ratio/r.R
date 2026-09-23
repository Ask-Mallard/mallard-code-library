# Fisher's exact test with the conditional maximum-likelihood odds ratio and its exact 95% interval
#
# For a 2x2 table with small expected counts. R's fisher.test reports the CONDITIONAL MLE of the
# odds ratio (conditioning on both margins), not the sample odds ratio ad/bc. The two differ, most
# visibly in small tables, and a report should say which it gives.

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d), all(d$exposed %in% c(0, 1)), all(d$outcome %in% c(0, 1)))

a <- sum(d$exposed == 1 & d$outcome == 1); b <- sum(d$exposed == 1 & d$outcome == 0)
c <- sum(d$exposed == 0 & d$outcome == 1); dd <- sum(d$exposed == 0 & d$outcome == 0)

# PINNED: the table is built explicitly, exposed in row 1 and the outcome in column 1, so the odds
# ratio is the odds of the outcome in the exposed over the unexposed. table(d$exposed, d$outcome)
# would sort both to 0 first; with an odd number of flips the odds ratio silently inverts.
tab <- matrix(c(a, b, c, dd), nrow = 2, byrow = TRUE,
              dimnames = list(exposed = c("1", "0"), outcome = c("1", "0")))
fit <- fisher.test(tab, conf.level = 0.95)

print(tab)
cat(sprintf("conditional MLE odds ratio %.3f, exact 95%% %.3f to %.3f, p = %.4g\n",
            fit$estimate, fit$conf.int[1], fit$conf.int[2], fit$p.value))
cat(sprintf("sample odds ratio ad/bc %.3f (not what fisher.test reports)\n", a * dd / (b * c)))

cat("\n--- HARNESS ---\n")
cat(sprintf("odds_ratio=%.10f\n", fit$estimate))
cat(sprintf("log_odds_ratio_lcl=%.10f\n", log(fit$conf.int[1])))
cat(sprintf("log_odds_ratio_ucl=%.10f\n", log(fit$conf.int[2])))
cat(sprintf("log_odds_ratio=%.10f\n", log(fit$estimate)))
cat(sprintf("sample_odds_ratio=%.10f\n", a * dd / (b * c)))
cat(sprintf("p_value=%.10e\n", fit$p.value))
cat(sprintf("n=%d\n", nrow(d)))
