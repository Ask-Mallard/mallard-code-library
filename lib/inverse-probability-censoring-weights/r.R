# Inverse probability of censoring weights (IPCW) for informative dropout.
#
# 1. DROPOUT THAT DEPENDS ON PROGNOSIS BIASES THE COMPLETERS. When sicker patients drop out more (and
#    more so in one arm), the patients whose outcome is observed are not representative of their arm.
#    The complete-case risk difference is biased even in a randomized trial.
#
# 2. THE WEIGHTS. A model for staying in the study given arm and the prognostic factor gives each
#    completer a weight 1 / P(stayed | A, L); weighted, the completers stand in for everyone randomized.
#    That is valid if dropout is missing at random given A and L (no unmeasured cause of both dropout
#    and outcome). Here the dropout model includes the arm-by-factor interaction as well.
#
# 3. THE SE COMES FROM THE STACKED EQUATIONS: the dropout model's score and the two weighted risks, so
#    the estimation of the weights is included. Treating the weights as known is reported beside it.
#
# 4. With time-to-event outcomes the same idea gives time-varying weights and a weighted Kaplan-Meier
#    or Cox model (see target-trial-clone-censor-weight).

d <- read.csv("fixture.csv")
d$stayed <- as.integer(!is.na(d$event))
cm <- glm(stayed ~ treated * l, data = d, family = binomial)
p <- fitted(cm)
X <- model.matrix(cm)
a <- d$treated
r <- d$stayed
y <- ifelse(r == 1, d$event, 0)
w <- r / p  # PINNED: 1 / P(stayed | A, L) for completers, 0 for dropouts

mu1 <- sum(w * a * y) / sum(w * a)
mu0 <- sum(w * (1 - a) * y) / sum(w * (1 - a))
rd <- mu1 - mu0

# Stacked equations: (dropout coefficients, mu1, mu0).
k <- ncol(X)
psi <- cbind(X * (r - p), w * a * (y - mu1), w * (1 - a) * (y - mu0))
A <- matrix(0, k + 2, k + 2)
A[1:k, 1:k] <- -crossprod(X * (p * (1 - p)), X)
# d/dgamma of r / p = -r (1 - p) / p x
A[k + 1, 1:k] <- colSums(-(r * (1 - p) / p) * a * (y - mu1) * X)
A[k + 1, k + 1] <- -sum(w * a)
A[k + 2, 1:k] <- colSums(-(r * (1 - p) / p) * (1 - a) * (y - mu0) * X)
A[k + 2, k + 2] <- -sum(w * (1 - a))
Ai <- solve(A)
V <- Ai %*% crossprod(psi) %*% t(Ai)
cvec <- c(rep(0, k), 1, -1)
se <- sqrt(drop(t(cvec) %*% V %*% cvec))
se_fixed <- sqrt(sum((w * a * (y - mu1))^2) / sum(w * a)^2 + sum((w * (1 - a) * (y - mu0))^2) / sum(w * (1 - a))^2)

cc <- mean(d$event[a == 1 & r == 1]) - mean(d$event[a == 0 & r == 1])

z <- qnorm(0.975)
cat(sprintf("IPCW risk difference %.4f (95%% CI %.4f to %.4f); complete-case %.4f\n", rd, rd - z * se, rd + z * se, cc))

cat("\n--- HARNESS ---\n")
cat(sprintf("risk_difference=%.10f\nrisk_difference_se=%.10f\n", rd, se))
cat(sprintf("risk_difference_lcl=%.10f\nrisk_difference_ucl=%.10f\n", rd - z * se, rd + z * se))
cat(sprintf("risk_difference_se_weights_fixed=%.10f\n", se_fixed))
cat(sprintf("complete_case_risk_difference=%.10f\nmax_weight=%.10f\n", cc, max(w)))
cat(sprintf("dropped=%d\nn=%d\n", sum(r == 0), nrow(d)))
