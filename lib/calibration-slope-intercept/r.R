# Calibration of a published prediction model in new patients: slope, intercept, O/E, and the curve.
#
# 1. DISCRIMINATION IS NOT CALIBRATION. The AUC is unchanged by any monotone distortion of the predicted
#    risks; a model can rank patients well and still predict 30% for patients whose risk is 15%.
#    Calibration asks whether predicted risks match observed ones.
#
# 2. THREE NUMBERS, EACH FROM ITS OWN MODEL (Van Calster et al. 2019):
#    - calibration SLOPE: the coefficient of the linear predictor in glm(y ~ lp). Below 1 means
#      predictions too extreme (overfitting); here the truth is 2/3.
#    - calibration-in-the-large (CITL): the intercept of glm(y ~ offset(lp)), with the slope FIXED at 1.
#      The intercept of the slope model is not CITL; it answers a different question. Below 0 means
#      over-prediction on average.
#    - observed/expected: total events over the sum of predicted risks.
#
# 3. THE CURVE: observed against predicted risk, smoothed (loess, or a restricted cubic spline in the
#    predicted logit), is the fuller picture and should be plotted. It is not a harness key: R's loess
#    and Python's lowess are different smoothers, so it cannot be compared across languages. The
#    Hosmer-Lemeshow test is not a substitute: its grouping is arbitrary and in large samples it rejects
#    trivial miscalibration.

d <- read.csv("fixture.csv")
slope_fit <- glm(event ~ lp, data = d, family = binomial)
citl_fit <- glm(event ~ offset(lp), data = d, family = binomial)
p <- plogis(d$lp)
oe <- sum(d$event) / sum(p)
auc <- {
  x <- d$lp[d$event == 1]; y <- d$lp[d$event == 0]
  mean(outer(x, y, ">") + 0.5 * outer(x, y, "=="))
}

# The SEs written out as inverse information (what vcov() gives), so both languages agree exactly.
X <- cbind(1, d$lp)
ps <- fitted(slope_fit)
slope_se <- sqrt(solve(crossprod(X * (ps * (1 - ps)), X))[2, 2])
pc <- plogis(coef(citl_fit)[[1]] + d$lp)
citl_se <- 1 / sqrt(sum(pc * (1 - pc)))

cat(sprintf("calibration slope %.3f, calibration-in-the-large %.3f, O/E %.3f, AUC %.3f\n",
            coef(slope_fit)[["lp"]], coef(citl_fit)[[1]], oe, auc))

cat("\n--- HARNESS ---\n")
cat(sprintf("calibration_slope=%.10f\ncalibration_slope_se=%.10f\n", coef(slope_fit)[["lp"]], slope_se))
cat(sprintf("calibration_in_the_large=%.10f\ncalibration_in_the_large_se=%.10f\n", coef(citl_fit)[[1]], citl_se))
cat(sprintf("slope_model_intercept=%.10f\n", coef(slope_fit)[[1]]))
cat(sprintf("observed_expected=%.10f\nauc=%.10f\n", oe, auc))
cat(sprintf("n=%d\nevents=%d\n", nrow(d), sum(d$event)))
