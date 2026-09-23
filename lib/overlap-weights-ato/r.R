# Overlap weights: the average treatment effect in the overlap population (ATO).
#
# 1. THE TARGET IS THE ATO, NOT THE ATE. Each treated patient is weighted by 1 - e and each control by
#    e, where e is the propensity score. Patients whose treatment was close to a coin toss count most;
#    patients almost certain to receive (or not receive) treatment count least. That is the
#    population with clinical equipoise, and when the effect varies it is a different number from
#    the ATE (here 2.25 against 2.00). Report it as the ATO.
#
# 2. NO EXTREME WEIGHTS. Overlap weights are bounded by 1, so there is nothing to trim or truncate.
#    That is their practical advantage over inverse-probability weights (Li, Morgan and Zaslavsky 2018).
#
# 3. EXACT BALANCE. With a logistic propensity model, overlap weights balance the MEAN of every
#    covariate in the model exactly: the weighted SMD of severity is zero to rounding. A non-zero
#    value means the weights or the model are not what this file says.
#
# 4. THE SE ACCOUNTS FOR THE ESTIMATED SCORE. The propensity model and the two weighted means are
#    solved as one stack of estimating equations, and the sandwich is taken over the whole stack.
#    Treating the weights as known (a robust SE from a weighted regression) is reported beside it;
#    for the ATO it is not guaranteed to be conservative.

d <- read.csv("fixture.csv")
a <- d$treated
y <- d$y
X <- cbind(1, d$severity, d$comorbid)

ps <- glm(treated ~ severity + comorbid, data = d, family = binomial)
e <- fitted(ps)
w <- ifelse(a == 1, 1 - e, e)  # PINNED: overlap weights

mu1 <- sum(w * a * y) / sum(w * a)
mu0 <- sum(w * (1 - a) * y) / sum(w * (1 - a))
ato <- mu1 - mu0

# PINNED: the stacked M-estimation sandwich over (propensity coefficients, mu1, mu0).
n <- nrow(d)
k <- ncol(X)
psi <- cbind(X * (a - e), w * a * (y - mu1), w * (1 - a) * (y - mu0))
de <- e * (1 - e)
A <- matrix(0, k + 2, k + 2)
A[1:k, 1:k] <- -crossprod(X * de, X)
A[k + 1, 1:k] <- colSums(-a * (y - mu1) * de * X)
A[k + 1, k + 1] <- -sum(w * a)
A[k + 2, 1:k] <- colSums((1 - a) * (y - mu0) * de * X)
A[k + 2, k + 2] <- -sum(w * (1 - a))
Ainv <- solve(A)
V <- Ainv %*% crossprod(psi) %*% t(Ainv)
contrast <- c(rep(0, k), 1, -1)
se <- sqrt(drop(t(contrast) %*% V %*% contrast))

# The weights treated as known: the HC0 sandwich from the weighted regression.
bread <- 1 / c(sum(w * a), sum(w * (1 - a)))
se_fixed <- sqrt(sum((w * a * (y - mu1))^2) * bread[1]^2 + sum((w * (1 - a) * (y - mu0))^2) * bread[2]^2)

wmean <- function(x, g) sum(w * g * x) / sum(w * g)
smd_w <- (wmean(d$severity, a) - wmean(d$severity, 1 - a)) /
  sqrt((var(d$severity[a == 1]) + var(d$severity[a == 0])) / 2)

z <- qnorm(0.975)
cat(sprintf("ATO %.3f (95%% CI %.3f to %.3f)\n", ato, ato - z * se, ato + z * se))

cat("\n--- HARNESS ---\n")
cat(sprintf("ato=%.10f\nato_se=%.10f\nato_lcl=%.10f\nato_ucl=%.10f\n", ato, se, ato - z * se, ato + z * se))
cat(sprintf("ato_se_weights_fixed=%.10f\n", se_fixed))
cat(sprintf("smd_severity_weighted=%.10f\n", smd_w))
cat(sprintf("max_weight=%.10f\nn=%d\n", max(w), n))
