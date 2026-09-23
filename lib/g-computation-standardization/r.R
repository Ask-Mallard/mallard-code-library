# G-computation (marginal standardization): a population risk difference and risk ratio from a
# logistic outcome model.
#
# 1. THE COEFFICIENT IS NOT THE ANSWER. exp(coef) is the CONDITIONAL odds ratio, the same in every
#    stratum of the covariates. The question a trial answers is marginal: what would the risk be if
#    everyone were treated, against no one? G-computation predicts both risks for every patient from
#    the fitted model and averages them.
#
# 2. THE ODDS RATIO IS NON-COLLAPSIBLE. The marginal odds ratio is smaller than the conditional one
#    even with no confounding at all (here log 0.61 against 0.70 in the data-generating process), so
#    reporting exp(coef) as "the effect of treatment" in the population overstates it. The risk
#    difference and risk ratio have no such problem once standardized.
#
# 3. THE SE COVERS BOTH SOURCES OF VARIATION. The model and the two standardized risks are solved as
#    one stack of estimating equations and the sandwich is taken over all of them (as R's stdReg
#    does). That includes the sampling of the covariates, which the population effect averages over.
#    The delta method on the model's covariance alone (the default of marginaleffects' avg_comparisons)
#    treats the covariates as fixed and is reported beside it; it is smaller (here only slightly,
#    because the model explains little of the outcome).

d <- read.csv("fixture.csv")
fit <- glm(event ~ treated + severity + comorbid, data = d, family = binomial)

d1 <- transform(d, treated = 1)
d0 <- transform(d, treated = 0)
p1 <- predict(fit, newdata = d1, type = "response")
p0 <- predict(fit, newdata = d0, type = "response")
mu1 <- mean(p1)
mu0 <- mean(p0)
rd <- mu1 - mu0

# PINNED: the stacked sandwich over (outcome-model coefficients, mu1, mu0).
X <- model.matrix(fit)
X1 <- model.matrix(fit, data = d1)
X0 <- model.matrix(fit, data = d0)
p <- fitted(fit)
n <- nrow(d)
k <- ncol(X)
psi <- cbind(X * (d$event - p), p1 - mu1, p0 - mu0)
A <- matrix(0, k + 2, k + 2)
A[1:k, 1:k] <- -crossprod(X * (p * (1 - p)), X)
A[k + 1, 1:k] <- colSums(X1 * (p1 * (1 - p1)))
A[k + 1, k + 1] <- -n
A[k + 2, 1:k] <- colSums(X0 * (p0 * (1 - p0)))
A[k + 2, k + 2] <- -n
Ainv <- solve(A)
V <- Ainv %*% crossprod(psi) %*% t(Ainv)
g_rd <- c(rep(0, k), 1, -1)
g_lrr <- c(rep(0, k), 1 / mu1, -1 / mu0)
rd_se <- sqrt(drop(t(g_rd) %*% V %*% g_rd))
lrr <- log(mu1 / mu0)
lrr_se <- sqrt(drop(t(g_lrr) %*% V %*% g_lrr))

# The delta method on vcov(fit) alone: covariates treated as fixed.
grad <- colMeans(X1 * (p1 * (1 - p1))) - colMeans(X0 * (p0 * (1 - p0)))
info_inv <- solve(crossprod(X * (p * (1 - p)), X))  # vcov(fit), written out so both languages agree exactly
rd_se_conditional <- sqrt(drop(t(grad) %*% info_inv %*% grad))

crude <- mean(d$event[d$treated == 1]) - mean(d$event[d$treated == 0])
z <- qnorm(0.975)
cat(sprintf("risk difference %.4f (95%% CI %.4f to %.4f); risk ratio %.3f\n", rd, rd - z * rd_se,
            rd + z * rd_se, exp(lrr)))

cat("\n--- HARNESS ---\n")
cat(sprintf("risk_treated=%.10f\nrisk_untreated=%.10f\n", mu1, mu0))
cat(sprintf("risk_difference=%.10f\nrisk_difference_se=%.10f\n", rd, rd_se))
cat(sprintf("risk_difference_lcl=%.10f\nrisk_difference_ucl=%.10f\n", rd - z * rd_se, rd + z * rd_se))
cat(sprintf("risk_difference_se_conditional=%.10f\n", rd_se_conditional))
cat(sprintf("log_risk_ratio=%.10f\nlog_risk_ratio_se=%.10f\n", lrr, lrr_se))
cat(sprintf("log_marginal_odds_ratio=%.10f\n", log(mu1 / (1 - mu1)) - log(mu0 / (1 - mu0))))
cat(sprintf("log_conditional_odds_ratio=%.10f\n", coef(fit)[["treated"]]))
cat(sprintf("crude_risk_difference=%.10f\nn=%d\n", crude, n))
