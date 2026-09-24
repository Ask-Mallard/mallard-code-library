# Random-effects meta-analysis of odds ratios: REML for tau^2, Hartung-Knapp for the interval.
#
# 1. RANDOM EFFECTS ESTIMATE THE MEAN OF A DISTRIBUTION OF TRUE EFFECTS, not one common effect. Report
#    tau^2 (how much the true effects vary) and a PREDICTION INTERVAL (where the effect in a new trial
#    would likely lie) beside the summary; I^2 is a proportion of variability, not an amount.
#
# 2. REML FOR tau^2 (not DerSimonian-Laird, which underestimates it with few trials), and HARTUNG-KNAPP
#    for the interval of the mean: a t distribution on k - 1 df with a rescaled variance. The normal
#    (z) interval with DerSimonian-Laird is too narrow when trials are few or heterogeneous; it is
#    reported beside it.
#
# 3. Pinned: the log odds ratio with Woolf variances (escalc measure = "OR"; no zero cells here, so no
#    continuity correction), REML to a convergence threshold of 1e-12, test = "knha".

suppressPackageStartupMessages(library(metafor))

d <- read.csv("fixture.csv")
es <- escalc(measure = "OR", ai = events_treated, bi = n_treated - events_treated,
             ci = events_control, di = n_control - events_control, data = d)
fit <- rma(yi, vi, data = es, method = "REML", test = "knha", control = list(threshold = 1e-12, maxiter = 1000))
pred <- predict(fit)
dl <- rma(yi, vi, data = es, method = "DL", test = "z")

cat(sprintf("summary OR %.3f (HK 95%% CI %.3f to %.3f); tau^2 %.4f, I^2 %.1f%%; prediction interval %.3f to %.3f\n",
            exp(coef(fit)), exp(fit$ci.lb), exp(fit$ci.ub), fit$tau2, fit$I2, exp(pred$pi.lb), exp(pred$pi.ub)))

cat("\n--- HARNESS ---\n")
cat(sprintf("log_or=%.10f\nlog_or_se=%.10f\nlog_or_lcl=%.10f\nlog_or_ucl=%.10f\n", coef(fit)[[1]], fit$se, fit$ci.lb, fit$ci.ub))
cat(sprintf("tau2=%.10f\ni2=%.10f\nq=%.10f\n", fit$tau2, fit$I2, fit$QE))
cat(sprintf("prediction_lcl=%.10f\nprediction_ucl=%.10f\n", pred$pi.lb, pred$pi.ub))
cat(sprintf("dl_log_or=%.10f\ndl_se_z=%.10f\ntrials=%d\n", coef(dl)[[1]], dl$se, fit$k))
