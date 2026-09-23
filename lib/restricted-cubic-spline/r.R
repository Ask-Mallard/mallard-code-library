# Restricted (natural) cubic spline for a nonlinear continuous exposure
#
# Fits a smooth curve for age instead of forcing a straight line or cutting age into categories. The
# individual spline coefficients have no meaning on their own and depend on the basis; report what the
# curve says: predicted values at chosen ages, and contrasts between them with intervals.

library(splines)

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d))

# PINNED: the knots are written out (30, 45, 60, 70, 80; the outer two are the boundary knots, beyond
# which the curve is linear). ns(age, df = 4) would place knots at the data's quantiles, and packages
# compute quantiles differently, so the same call can fit different curves in different software.
basis <- function(age) ns(age, knots = c(45, 60, 70), Boundary.knots = c(30, 80))
fit <- lm(sbp ~ basis(age), data = d)

new <- data.frame(age = c(40, 50, 70))
Xn <- model.matrix(~ basis(age), new)
pred <- as.vector(Xn %*% coef(fit))
contrast <- Xn[3, ] - Xn[2, ]              # age 70 minus age 50
est <- sum(contrast * coef(fit))
se <- sqrt(as.numeric(t(contrast) %*% vcov(fit) %*% contrast))
tq <- qt(0.975, df.residual(fit))

cat(sprintf("predicted SBP at 40, 50, 70: %.2f, %.2f, %.2f\n", pred[1], pred[2], pred[3]))
cat(sprintf("age 70 vs 50: %.2f mmHg, 95%% %.2f to %.2f\n", est, est - tq * se, est + tq * se))

cat("\n--- HARNESS ---\n")
cat(sprintf("pred_40=%.10f\npred_50=%.10f\npred_70=%.10f\n", pred[1], pred[2], pred[3]))
cat(sprintf("contrast_70_50=%.10f\ncontrast_70_50_se=%.10f\n", est, se))
cat(sprintf("contrast_70_50_lcl=%.10f\ncontrast_70_50_ucl=%.10f\n", est - tq * se, est + tq * se))
cat(sprintf("n=%d\n", nrow(d)))
