# Sensitivity to unmeasured confounding: the E-value and a simple quantitative bias analysis.
#
# 1. THE E-VALUE (VanderWeele and Ding 2017) is the minimum strength of association, on the risk-ratio
#    scale, that an unmeasured confounder would need with BOTH the exposure and the outcome to explain
#    away the observed risk ratio: RR + sqrt(RR (RR - 1)) for RR > 1. Report it for the point estimate
#    and for the confidence limit closer to 1. It is a summary of robustness, not a bias correction, and
#    it says nothing about whether such a confounder exists.
#
# 2. A QUANTITATIVE BIAS ANALYSIS goes further when the confounder can be named and its associations
#    estimated from external evidence: the prevalence of the confounder among the exposed (p1) and
#    unexposed (p0), and its risk ratio with the outcome (RR_UY). The classic bias factor
#    (p1 (RR_UY - 1) + 1) / (p0 (RR_UY - 1) + 1) divides the observed RR. It assumes no interaction
#    between the confounder and the exposure. The values below are this fixture's true ones, so the
#    adjusted RR should recover the causal 1.5; in practice they come from the literature, and a range
#    of values (or a probabilistic analysis, episensr) is reported.
#
# 3. The observed RR here is the crude 2x2 risk ratio with a log-scale (Katz) interval.

library(EValue)

d <- read.csv("fixture.csv")
a <- sum(d$event[d$exposed == 1]); n1 <- sum(d$exposed == 1)
c0 <- sum(d$event[d$exposed == 0]); n0 <- sum(d$exposed == 0)
log_rr <- log((a / n1) / (c0 / n0))
se <- sqrt(1 / a - 1 / n1 + 1 / c0 - 1 / n0)
z <- qnorm(0.975)
rr <- exp(log_rr); lo <- exp(log_rr - z * se); hi <- exp(log_rr + z * se)

ev <- evalues.RR(est = rr, lo = lo, hi = hi)

# PINNED: the bias parameters (external evidence in practice; the fixture's true values here).
p1 <- 0.5236
p0 <- 0.2879
rr_uy <- 2.5
bias <- (p1 * (rr_uy - 1) + 1) / (p0 * (rr_uy - 1) + 1)

cat(sprintf("observed RR %.2f (95%% CI %.2f to %.2f); E-value %.2f (CI limit %.2f); bias-adjusted RR %.2f\n",
            rr, lo, hi, ev["E-values", "point"], ev["E-values", "lower"], rr / bias))

cat("\n--- HARNESS ---\n")
cat(sprintf("log_rr_observed=%.10f\nlog_rr_observed_se=%.10f\n", log_rr, se))
cat(sprintf("e_value_point=%.10f\ne_value_ci=%.10f\n", ev["E-values", "point"], ev["E-values", "lower"]))
cat(sprintf("bias_factor=%.10f\nlog_rr_bias_adjusted=%.10f\n", bias, log_rr - log(bias)))
cat(sprintf("n=%d\n", nrow(d)))
