# Delta-adjusted tipping-point analysis for missing outcomes (a pattern-mixture sensitivity analysis).
#
# 1. THE PRIMARY ANALYSIS ASSUMES MISSING AT RANDOM, and that cannot be checked. Here each arm's outcome
#    is modelled on the baseline covariate among patients with an outcome, the model predicts the missing
#    outcomes, and each arm's mean averages observed and predicted values (regression imputation of the
#    mean). Under MAR that is unbiased; the complete-case difference is not.
#
# 2. THE SENSITIVITY ANALYSIS SHIFTS THE MISSING. A shift delta is added to the predicted outcomes of
#    the TREATED patients whose outcome is missing: delta < 0 says they did worse than MAR predicts
#    (missing not at random). The TIPPING POINT is the delta at which the lower 95% limit of the effect
#    reaches zero. Clinicians then judge whether a shift that large is plausible. Shifting both arms, or
#    a grid of (delta_treated, delta_control), are common variants.
#
# 3. THE STANDARD ERROR at each delta is from the stacked estimating equations (each arm's outcome model
#    and its mean), so the estimation of the imputation model is included. Written out, so both
#    languages give the same number; the usual alternative is multiple imputation with a delta added to
#    the imputed values (for example mice's post-processing), which adds Monte Carlo error.

d <- read.csv("fixture.csv")
d$r <- as.integer(!is.na(d$y))

arm_fit <- function(a) {
  s <- d[d$treated == a, ]
  X <- cbind(1, s$x)
  obs <- s$r == 1
  b <- solve(crossprod(X[obs, ]), crossprod(X[obs, ], s$y[obs]))
  list(s = s, X = X, obs = obs, b = drop(b))
}
arms <- list(`1` = arm_fit(1), `0` = arm_fit(0))

effect_at <- function(delta) {
  est <- c()
  var <- 0
  for (a in c("1", "0")) {
    f <- arms[[a]]
    dl <- if (a == "1") delta else 0
    pred <- drop(f$X %*% f$b) + dl
    val <- ifelse(f$obs, f$s$y, pred)
    mu <- mean(val)
    n <- nrow(f$X)
    # Stacked equations (b0, b1, mu) for this arm; the arms are independent.
    res <- ifelse(f$obs, f$s$y - drop(f$X %*% f$b), 0)
    psi <- cbind(f$X * res, val - mu)
    A <- matrix(0, 3, 3)
    A[1:2, 1:2] <- -crossprod(f$X[f$obs, ])
    A[3, 1:2] <- colSums(f$X * (!f$obs))
    A[3, 3] <- -n
    Ai <- solve(A)
    V <- Ai %*% crossprod(psi) %*% t(Ai)
    est[a] <- mu
    var <- var + V[3, 3]
  }
  c(effect = est[["1"]] - est[["0"]], se = sqrt(var))
}

z <- qnorm(0.975)
mar <- effect_at(0)
tip <- uniroot(function(dl) { e <- effect_at(dl); e[["effect"]] - z * e[["se"]] }, c(-50, 0), tol = 1e-12)$root

cc <- mean(d$y[d$treated == 1], na.rm = TRUE) - mean(d$y[d$treated == 0], na.rm = TRUE)

cat(sprintf("effect under MAR %.3f (95%% CI %.3f to %.3f); tipping point: treated missing outcomes %.2f lower\n",
            mar[["effect"]], mar[["effect"]] - z * mar[["se"]], mar[["effect"]] + z * mar[["se"]], -tip))

cat("\n--- HARNESS ---\n")
cat(sprintf("effect=%.10f\neffect_se=%.10f\n", mar[["effect"]], mar[["se"]]))
cat(sprintf("effect_lcl=%.10f\neffect_ucl=%.10f\n", mar[["effect"]] - z * mar[["se"]], mar[["effect"]] + z * mar[["se"]]))
cat(sprintf("tipping_delta=%.10f\n", tip))
cat(sprintf("effect_at_delta_minus2=%.10f\n", effect_at(-2)[["effect"]]))
cat(sprintf("complete_case_difference=%.10f\nmissing=%d\nn=%d\n", cc, sum(d$r == 0), nrow(d)))
