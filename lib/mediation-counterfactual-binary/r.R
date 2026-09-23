# Causal mediation with a binary outcome: natural direct and indirect effects by the mediation formula.
#
# 1. BARON-KENNY AND PRODUCT-OF-COEFFICIENTS ARE INVALID HERE. With a logistic outcome model the odds
#    ratio is non-collapsible, so the change in the treatment coefficient when the mediator is added is
#    not an indirect effect, and a product of coefficients has no counterfactual meaning. Both also
#    ignore the treatment-by-mediator interaction. The difference-method number is reported beside the
#    counterfactual one so the gap is visible.
#
# 2. THE COUNTERFACTUAL DEFINITIONS (Robins-Greenland, Pearl). With E[Y(a, M(a*))] the risk under
#    treatment a with the mediator as it would be under a*:
#      natural direct effect   NDE = E[Y(1, M(0))] - E[Y(0, M(0))]
#      natural indirect effect NIE = E[Y(1, M(1))] - E[Y(1, M(0))]
#    and NDE + NIE = the total effect. The mediation formula computes each from a mediator model and an
#    outcome model WITH the interaction, averaged over the covariates (VanderWeele 2015). With a binary
#    mediator the sum over its two values is exact; no simulation is needed, so the result is
#    reproducible. R's mediation package computes the same quantities by quasi-Bayesian simulation.
#
# 3. THE ASSUMPTIONS ARE STRONGER THAN FOR A TOTAL EFFECT: no unmeasured confounding of
#    treatment-outcome, mediator-outcome and treatment-mediator, and no mediator-outcome confounder that
#    is itself affected by treatment. A sensitivity analysis for the mediator-outcome assumption belongs
#    in the plan.
#
# 4. THE SE COMES FROM ONE STACK: both models' scores and the three potential risks, with the sandwich
#    over all of them. Power for mediation needs its own simulation.

d <- read.csv("fixture.csv")
n <- nrow(d)
mfit <- glm(mediator ~ treated + c, data = d, family = binomial)
yfit <- glm(event ~ treated * mediator + c, data = d, family = binomial)
g <- coef(mfit)
b <- coef(yfit)

XM <- function(a) cbind(1, a, d$c)
XY <- function(a, m) cbind(1, a, m, d$c, a * m)  # column order of coef(yfit)
stopifnot(identical(names(b), c("(Intercept)", "treated", "mediator", "c", "treated:mediator")))
expit <- plogis
pm <- function(a) drop(expit(XM(a) %*% g))
py <- function(a, m) drop(expit(XY(a, m) %*% b))
mu_i <- function(a, as) py(a, 1) * pm(as) + py(a, 0) * (1 - pm(as))

# PINNED: the mediation formula, summed exactly over the binary mediator.
mu <- c(mu10 = mean(mu_i(1, 0)), mu00 = mean(mu_i(0, 0)), mu11 = mean(mu_i(1, 1)))
nde <- mu[["mu10"]] - mu[["mu00"]]
nie <- mu[["mu11"]] - mu[["mu10"]]

# PINNED: the stacked sandwich over (mediator coefficients, outcome coefficients, three risks).
Xm <- model.matrix(mfit)
Xy <- model.matrix(yfit)
p_m <- fitted(mfit)
p_y <- fitted(yfit)
k1 <- ncol(Xm)
k2 <- ncol(Xy)
pairs <- list(c(1, 0), c(0, 0), c(1, 1))
psi <- cbind(Xm * (d$mediator - p_m), Xy * (d$event - p_y),
             sapply(seq_along(pairs), function(j) mu_i(pairs[[j]][1], pairs[[j]][2]) - mu[j]))
K <- k1 + k2 + 3
A <- matrix(0, K, K)
A[1:k1, 1:k1] <- -crossprod(Xm * (p_m * (1 - p_m)), Xm)
A[k1 + 1:k2, k1 + 1:k2] <- -crossprod(Xy * (p_y * (1 - p_y)), Xy)
for (j in seq_along(pairs)) {
  a <- pairs[[j]][1]; as <- pairs[[j]][2]
  q <- pm(as); y1 <- py(a, 1); y0 <- py(a, 0)
  A[k1 + k2 + j, 1:k1] <- colSums(XM(as) * ((y1 - y0) * q * (1 - q)))
  A[k1 + k2 + j, k1 + 1:k2] <- colSums(XY(a, 1) * (y1 * (1 - y1) * q) + XY(a, 0) * (y0 * (1 - y0) * (1 - q)))
  A[k1 + k2 + j, k1 + k2 + j] <- -n
}
Ainv <- solve(A)
V <- Ainv %*% crossprod(psi) %*% t(Ainv)
cvec <- function(w) c(rep(0, k1 + k2), w)
se <- function(w) sqrt(drop(t(cvec(w)) %*% V %*% cvec(w)))
nde_se <- se(c(1, -1, 0))
nie_se <- se(c(-1, 0, 1))

# The difference method on the log-odds scale, for comparison only: not an indirect effect here.
total_coef <- coef(glm(event ~ treated + c, data = d, family = binomial))[["treated"]]
direct_coef <- coef(glm(event ~ treated + mediator + c, data = d, family = binomial))[["treated"]]
nie_log_or <- qlogis(mu[["mu11"]]) - qlogis(mu[["mu10"]])

z <- qnorm(0.975)
cat(sprintf("NDE %.4f (95%% CI %.4f to %.4f); NIE %.4f (%.4f to %.4f); proportion mediated %.2f\n",
            nde, nde - z * nde_se, nde + z * nde_se, nie, nie - z * nie_se, nie + z * nie_se, nie / (nde + nie)))

cat("\n--- HARNESS ---\n")
cat(sprintf("natural_direct_effect=%.10f\nnatural_direct_effect_se=%.10f\n", nde, nde_se))
cat(sprintf("natural_indirect_effect=%.10f\nnatural_indirect_effect_se=%.10f\n", nie, nie_se))
cat(sprintf("total_effect=%.10f\ntotal_effect_se=%.10f\n", nde + nie, se(c(0, -1, 1))))
cat(sprintf("nie_log_odds_ratio=%.10f\ndifference_method_log_odds_ratio=%.10f\n", nie_log_or, total_coef - direct_coef))
cat(sprintf("n=%d\n", n))
