# Agreement between two raters on an ordinal scale: Cohen's kappa and quadratic-weighted kappa.
#
# 1. PERCENT AGREEMENT IS NOT AGREEMENT BEYOND CHANCE. Two raters who both call most patients "grade 1"
#    agree often by chance alone. Kappa = (observed - chance) / (1 - chance).
#
# 2. ORDINAL SCALES NEED WEIGHTS. Unweighted kappa counts a 1-versus-2 disagreement the same as 1-versus-4.
#    Quadratic weights, 1 - (i - j)^2 / (k - 1)^2, give partial credit; weighted kappa with quadratic
#    weights is close to an intraclass correlation. Pinned: both reported, the weights named (irr's
#    weight = "squared").
#
# 3. THE STANDARD ERROR used for an interval is the large-sample one of Fleiss, Cohen and Everitt (1969),
#    written out here; irr::kappa2 reports only a test of kappa = 0, whose SE is a different quantity.
#
# 4. Kappa depends on the prevalence of each category; report the marginal distributions beside it.

suppressPackageStartupMessages(library(irr))

d <- read.csv("fixture.csv")
lev <- 1:4
k <- length(lev)
n <- nrow(d)
p <- table(factor(d$rater_a, lev), factor(d$rater_b, lev)) / n
pr <- rowSums(p); pc <- colSums(p)

fleiss_se <- function(w) {
  po <- sum(w * p); pe <- sum(w * outer(pr, pc)); kap <- (po - pe) / (1 - pe)
  wr <- drop(w %*% pc); wc <- drop(t(w) %*% pr)
  inner <- sum(p * (w - outer(wr, wc, "+") * (1 - kap))^2)
  sqrt((inner - (kap - pe * (1 - kap))^2) / (n * (1 - pe)^2))
}
w0 <- diag(k)
w2 <- 1 - outer(lev, lev, "-")^2 / (k - 1)^2

k0 <- kappa2(d[, c("rater_a", "rater_b")], weight = "unweighted")$value
k2 <- kappa2(d[, c("rater_a", "rater_b")], weight = "squared")$value
se0 <- fleiss_se(w0); se2 <- fleiss_se(w2)

z <- qnorm(0.975)
cat(sprintf("percent agreement %.3f; kappa %.3f (%.3f to %.3f); quadratic-weighted kappa %.3f (%.3f to %.3f)\n",
            sum(diag(p)), k0, k0 - z * se0, k0 + z * se0, k2, k2 - z * se2, k2 + z * se2))

cat("\n--- HARNESS ---\n")
cat(sprintf("percent_agreement=%.10f\n", sum(diag(p))))
cat(sprintf("kappa=%.10f\nkappa_se=%.10f\n", k0, se0))
cat(sprintf("kappa_quadratic=%.10f\nkappa_quadratic_se=%.10f\n", k2, se2))
cat(sprintf("n=%d\n", n))
