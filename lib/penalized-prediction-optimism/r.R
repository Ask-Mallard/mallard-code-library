# Developing a penalized prediction model and correcting its apparent performance by the bootstrap.
#
# 1. APPARENT PERFORMANCE IS OPTIMISTIC. A model fitted and assessed on the same patients has a higher
#    AUC and a calibration slope closer to 1 than it will have in new patients.
#
# 2. SHRINK AT DEVELOPMENT. Ridge (L2-penalized) logistic regression on standardized predictors, the
#    intercept unpenalized, the penalty chosen by 5-fold cross-validated deviance over a fixed grid. The
#    penalty here is on the log-likelihood scale: maximize loglik - (lambda / 2) sum(beta^2). glmnet
#    optimizes the same kind of objective with its own scaling of lambda (divided by n) and internal
#    standardization; this entry writes the fit out (Newton-Raphson to 1e-12) so both languages agree.
#
# 3. CORRECT BY THE BOOTSTRAP, REPEATING EVERY STEP (Harrell). For each resample: choose the penalty by
#    cross-validation inside the resample, fit, and record the model's AUC and calibration slope in the
#    resample (apparent) and in the original data (test). The mean difference is the optimism; the
#    corrected estimate is the apparent one minus it. The 50 resamples and the folds come from the
#    fixture so both languages use the same ones; in practice use 200 or more.
#
# 4. The fixture includes an external validation cohort, so the corrected estimates can be compared with
#    what actually happens in new patients.

d <- read.csv("fixture.csv")
dev <- d[d$cohort == "development", ]
val <- d[d$cohort == "validation", ]
xs <- paste0("x", 1:10)
boot <- lapply(strsplit(read.csv("bootstrap.csv")$rows, " "), function(v) as.integer(v) + 1L)
lambdas <- exp(seq(log(0.1), log(100), length.out = 12))

ridge <- function(X, y, lambda) {
  mu <- colMeans(X); s <- apply(X, 2, sd)
  Z <- cbind(1, sweep(sweep(X, 2, mu), 2, s, "/"))
  pen <- diag(c(0, rep(lambda, ncol(X))))
  b <- rep(0, ncol(Z))
  repeat {
    p <- plogis(drop(Z %*% b))
    step <- solve(crossprod(Z * (p * (1 - p)), Z) + pen, crossprod(Z, y - p) - pen %*% b)
    b <- b + drop(step)
    if (max(abs(step)) < 1e-12) break
  }
  beta <- b[-1] / s
  c(b[1] - sum(beta * mu), beta)  # intercept and slopes on the original scale
}
lp <- function(coef, X) drop(cbind(1, X) %*% coef)
auc <- function(eta, y) { x <- eta[y == 1]; z <- eta[y == 0]; mean(outer(x, z, ">") + 0.5 * outer(x, z, "==")) }
slope <- function(eta, y) coef(glm(y ~ eta, family = binomial, control = glm.control(epsilon = 1e-14, maxit = 100)))[[2]]
deviance <- function(eta, y) -2 * sum(y * eta - log1p(exp(eta)))

choose_lambda <- function(X, y, fold) {
  cv <- sapply(lambdas, function(l) sum(sapply(sort(unique(fold)), function(f) {
    tr <- fold != f
    deviance(lp(ridge(X[tr, , drop = FALSE], y[tr], l), X[!tr, , drop = FALSE]), y[!tr])
  })))
  lambdas[which.min(cv)]
}

X <- as.matrix(dev[, xs]); y <- dev$event
lam <- choose_lambda(X, y, dev$fold)
fit <- ridge(X, y, lam)
app_auc <- auc(lp(fit, X), y); app_slope <- slope(lp(fit, X), y)

opt <- t(sapply(boot, function(idx) {
  Xb <- X[idx, , drop = FALSE]; yb <- y[idx]
  fb <- ridge(Xb, yb, choose_lambda(Xb, yb, ((seq_along(idx) - 1) %% 5) + 1))
  c(auc(lp(fb, Xb), yb) - auc(lp(fb, X), y), slope(lp(fb, Xb), yb) - slope(lp(fb, X), y))
}))
cor_auc <- app_auc - mean(opt[, 1]); cor_slope <- app_slope - mean(opt[, 2])
Xv <- as.matrix(val[, xs])
ext_auc <- auc(lp(fit, Xv), val$event); ext_slope <- slope(lp(fit, Xv), val$event)

cat(sprintf("lambda %.3f; AUC apparent %.3f, corrected %.3f, external %.3f; slope apparent %.3f, corrected %.3f, external %.3f\n",
            lam, app_auc, cor_auc, ext_auc, app_slope, cor_slope, ext_slope))

cat("\n--- HARNESS ---\n")
cat(sprintf("lambda=%.10f\n", lam))
cat(sprintf("apparent_auc=%.10f\noptimism_auc=%.10f\ncorrected_auc=%.10f\nexternal_auc=%.10f\n", app_auc, mean(opt[, 1]), cor_auc, ext_auc))
cat(sprintf("apparent_slope=%.10f\noptimism_slope=%.10f\ncorrected_slope=%.10f\nexternal_slope=%.10f\n", app_slope, mean(opt[, 2]), cor_slope, ext_slope))
cat(sprintf("corrected_minus_external_auc=%.10f\ncorrected_minus_external_slope=%.10f\n", cor_auc - ext_auc, cor_slope - ext_slope))
cat(sprintf("coef_x1=%.10f\nn_development=%d\n", fit[2], nrow(dev)))
