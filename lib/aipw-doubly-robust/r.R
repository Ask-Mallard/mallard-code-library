# Augmented inverse-probability weighting (AIPW): a doubly robust average treatment effect.
#
# 1. TWO MODELS, ONE NEEDS TO BE RIGHT. AIPW combines a propensity model and an outcome model. The
#    estimate is consistent if EITHER is correctly specified, which is the "doubly robust" property.
#    Inverse-probability weighting alone needs the propensity model right; g-computation alone needs
#    the outcome model right. Both are reported beside AIPW.
#
# 2. THE FORMULA, WRITTEN OUT. For each patient, the outcome model's predictions m1 and m0 are
#    corrected by the weighted residual of the arm they were actually in:
#      mu1 = mean( m1 + A (Y - m1) / e ),   mu0 = mean( m0 + (1 - A) (Y - m0) / (1 - e) ).
#    Unstabilized weights, no trimming. Packages (AIPW, tmle, WeightIt) add cross-fitting and
#    machine-learning nuisance models; the estimator underneath is this one.
#
# 3. THE SE IS FROM THE INFLUENCE FUNCTION: sd(IF) / sqrt(n), with
#    IF = m1 + A (Y - m1) / e - m0 - (1 - A) (Y - m0) / (1 - e) - ATE.
#    It ignores that the two models were estimated, which is valid when both are correctly specified
#    (AIPW is then efficient); when only one is right it can be off, and a bootstrap is the check.

d <- read.csv("fixture.csv")
a <- d$treated
y <- d$event

ps <- glm(treated ~ severity + comorbid, data = d, family = binomial)
e <- fitted(ps)
om <- glm(event ~ treated + severity + comorbid, data = d, family = binomial)
m1 <- predict(om, newdata = transform(d, treated = 1), type = "response")
m0 <- predict(om, newdata = transform(d, treated = 0), type = "response")

# PINNED: the AIPW estimating equation, unstabilized, untrimmed.
phi1 <- m1 + a * (y - m1) / e
phi0 <- m0 + (1 - a) * (y - m0) / (1 - e)
ate <- mean(phi1) - mean(phi0)
n <- nrow(d)
IF <- phi1 - phi0 - ate
se <- sqrt(sum(IF^2)) / n

ipw <- mean(a * y / e) - mean((1 - a) * y / (1 - e))  # Horvitz-Thompson, for comparison
gcomp <- mean(m1) - mean(m0)

z <- qnorm(0.975)
cat(sprintf("AIPW risk difference %.4f (95%% CI %.4f to %.4f)\n", ate, ate - z * se, ate + z * se))

cat("\n--- HARNESS ---\n")
cat(sprintf("ate=%.10f\nate_se=%.10f\nate_lcl=%.10f\nate_ucl=%.10f\n", ate, se, ate - z * se, ate + z * se))
cat(sprintf("ipw_only=%.10f\ngcomp_only=%.10f\n", ipw, gcomp))
cat(sprintf("n=%d\n", n))
