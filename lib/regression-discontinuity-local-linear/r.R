# Sharp regression discontinuity: local linear estimation with robust bias-corrected inference.
#
# 1. THE EFFECT IS THE JUMP AT THE CUT-OFF, for patients near it. It says nothing directly about
#    patients far from the threshold.
#
# 2. LOCAL, NOT GLOBAL. A straight line (or a high-order polynomial) fitted to the whole range on each
#    side misestimates the jump when the outcome curves (Gelman and Imbens 2019). Local linear fits
#    within a bandwidth h of the cut-off, weighted by a triangular kernel, are the standard.
#
# 3. THE BANDWIDTH IS CHOSEN BY A RULE, AND THE INTERVAL IS BIAS-CORRECTED. rdrobust picks h to
#    minimise the mean squared error of the point estimate ("mserd"). At that h the conventional
#    interval is too narrow, because the estimate still carries smoothing bias; the ROBUST
#    bias-corrected interval (Calonico, Cattaneo and Titiunik 2014) removes the bias with a local
#    quadratic fit over a second bandwidth b and widens the interval for it. Report the conventional
#    estimate with the robust interval, as rdrobust prints them.
#
# 4. EVERY OPTION NAMED: p = 1 (local linear), q = 2 (the bias correction), triangular kernel, one
#    common MSE-optimal bandwidth, nearest-neighbour variance with 3 matches.

library(rdrobust)

d <- read.csv("fixture.csv")
r <- rdrobust(d$y, d$score, c = 0, p = 1, q = 2, kernel = "triangular", bwselect = "mserd",
              vce = "nn", nnmatch = 3)

est <- r$coef["Conventional", 1]
lcl <- r$ci["Robust", "CI Lower"]
ucl <- r$ci["Robust", "CI Upper"]

# A global straight line on each side, for comparison: what the curvature does to it.
global <- coef(lm(y ~ treated * score, data = d))[["treated"]]

cat(sprintf("jump at the cut-off %.3f; robust 95%% CI %.3f to %.3f; bandwidth %.3f\n", est, lcl, ucl, r$bws["h", "left"]))

cat("\n--- HARNESS ---\n")
cat(sprintf("jump=%.10f\njump_se_conventional=%.10f\n", est, r$se["Conventional", 1]))
cat(sprintf("jump_bias_corrected=%.10f\njump_se_robust=%.10f\n", r$coef["Bias-Corrected", 1], r$se["Robust", 1]))
cat(sprintf("robust_lcl=%.10f\nrobust_ucl=%.10f\n", lcl, ucl))
cat(sprintf("conventional_lcl=%.10f\nconventional_ucl=%.10f\n", r$ci["Conventional", "CI Lower"], r$ci["Conventional", "CI Upper"]))
cat(sprintf("bandwidth_h=%.10f\nbandwidth_b=%.10f\n", r$bws["h", "left"], r$bws["b", "left"]))
cat(sprintf("n_left_in_h=%d\nn_right_in_h=%d\n", r$N_h[1], r$N_h[2]))
cat(sprintf("global_linear_jump=%.10f\nn=%d\n", global, nrow(d)))
