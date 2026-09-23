# Pearson and Spearman correlation, each with a 95% interval
#
# Pearson measures LINEAR association; Spearman measures monotonic association on ranks and resists
# outliers. Neither is a measure of agreement between two methods (see Bland-Altman), and neither says
# one variable causes the other.

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d))
n <- nrow(d)
z <- qnorm(0.975)

# Pearson: cor.test's interval uses Fisher's z transformation.
pear <- cor.test(d$bmi, d$sbp, method = "pearson", conf.level = 0.95)

# Spearman: cor.test gives NO interval for Spearman, so it is built on Fisher's z with the Fieller,
# Hartley and Pearson (1957) standard error. PINNED: 1.06 / sqrt(n - 3), not Pearson's 1 / sqrt(n - 3).
rho <- cor(d$bmi, d$sbp, method = "spearman")
se_s <- 1.06 / sqrt(n - 3)
sp_ci <- tanh(atanh(rho) + c(-1, 1) * z * se_s)

cat(sprintf("Pearson r %.4f, 95%% %.4f to %.4f\n", pear$estimate, pear$conf.int[1], pear$conf.int[2]))
cat(sprintf("Spearman rho %.4f, 95%% %.4f to %.4f\n", rho, sp_ci[1], sp_ci[2]))

cat("\n--- HARNESS ---\n")
cat(sprintf("pearson_r=%.10f\npearson_lcl=%.10f\npearson_ucl=%.10f\n", pear$estimate, pear$conf.int[1], pear$conf.int[2]))
cat(sprintf("spearman_rho=%.10f\nspearman_lcl=%.10f\nspearman_ucl=%.10f\n", rho, sp_ci[1], sp_ci[2]))
cat(sprintf("n=%d\n", n))
