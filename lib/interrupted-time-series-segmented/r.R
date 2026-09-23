# Interrupted time series: segmented regression with autocorrelation-robust (Newey-West) errors.
#
# 1. TWO EFFECTS, NOT ONE. The model has a pre-intervention trend, a LEVEL change at the start, and a
#    SLOPE change after it: rate ~ month + post + months_since. Report both; a level change alone can
#    hide a slope that erodes it.
#
# 2. NEIGHBOURING MONTHS ARE CORRELATED. Ordinary least squares treats them as independent and its
#    standard errors are too small when the correlation is positive, as it usually is. The Newey-West
#    (HAC) covariance allows correlation up to a stated lag. Pinned: lag 3, Bartlett weights, no
#    pre-whitening and no small-sample adjustment (NeweyWest's defaults choose the lag from the data
#    and pre-whiten, which statsmodels does not do). The lag is a choice to state in the protocol.
#    Prais-Winsten or an ARIMA model are alternatives that model the correlation instead.
#
# 3. SEASONALITY AND THE NUMBER OF POINTS. A monthly series with a seasonal cycle needs seasonal terms;
#    this fixture has none. Under about 8 points on either side the trends are not estimable with any
#    confidence.

library(sandwich)

d <- read.csv("fixture.csv")
fit <- lm(rate ~ month + post + months_since, data = d)

V <- NeweyWest(fit, lag = 3, prewhite = FALSE, adjust = FALSE)
se <- sqrt(diag(V))
se_ols <- sqrt(diag(vcov(fit)))
r <- resid(fit)
rho1 <- sum(r[-1] * r[-length(r)]) / sum(r^2)

z <- qnorm(0.975)
b <- coef(fit)
cat(sprintf("level change %.2f (95%% CI %.2f to %.2f); slope change %.3f a month (%.3f to %.3f)\n",
            b[["post"]], b[["post"]] - z * se[["post"]], b[["post"]] + z * se[["post"]],
            b[["months_since"]], b[["months_since"]] - z * se[["months_since"]], b[["months_since"]] + z * se[["months_since"]]))

cat("\n--- HARNESS ---\n")
cat(sprintf("level_change=%.10f\nlevel_change_se=%.10f\n", b[["post"]], se[["post"]]))
cat(sprintf("slope_change=%.10f\nslope_change_se=%.10f\n", b[["months_since"]], se[["months_since"]]))
cat(sprintf("pre_slope=%.10f\n", b[["month"]]))
cat(sprintf("level_change_se_ols=%.10f\nslope_change_se_ols=%.10f\n", se_ols[["post"]], se_ols[["months_since"]]))
cat(sprintf("residual_lag1_autocorrelation=%.10f\nn=%d\n", rho1, nrow(d)))
