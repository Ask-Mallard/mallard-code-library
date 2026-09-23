# Controlled interrupted time series: an intervention site against a control site.
#
# 1. WHY A CONTROL SERIES. A single-series ITS attributes everything that changes at the start date to
#    the intervention. Anything else that changes then (a coding change, a regional campaign, a new
#    guideline) is absorbed into the estimate. A control series exposed to the same other changes but not
#    to the intervention removes them. Here a drop of 2 hits both sites at the same month; the
#    single-series estimate is about -6 and the controlled one about -4.
#
# 2. THE MODEL, AS THE DIFFERENCE SERIES. With both sites observed on the same months, segmented
#    regression of the DIFFERENCE (intervention minus control) gives exactly the group-by-period
#    coefficients of the full interaction model on the stacked data (the control file checks it), and
#    its residuals are one ordered series, so the Newey-West covariance applies directly. Pinned: lag 3,
#    Bartlett, no pre-whitening, no adjustment.
#
# 3. THE ASSUMPTION: without the intervention, the difference between the sites would have continued its
#    pre-period trend. Choose a control that shares the intervention site's other exposures.

library(sandwich)

d <- read.csv("fixture.csv")
d$difference <- d$intervention_rate - d$control_rate
fit <- lm(difference ~ month + post + months_since, data = d)
se <- sqrt(diag(NeweyWest(fit, lag = 3, prewhite = FALSE, adjust = FALSE)))
b <- coef(fit)

# The single-series analysis of the intervention site, for comparison.
single <- coef(lm(intervention_rate ~ month + post + months_since, data = d))

z <- qnorm(0.975)
cat(sprintf("level change relative to control %.2f (95%% CI %.2f to %.2f); slope change %.3f\n",
            b[["post"]], b[["post"]] - z * se[["post"]], b[["post"]] + z * se[["post"]], b[["months_since"]]))

cat("\n--- HARNESS ---\n")
cat(sprintf("level_change=%.10f\nlevel_change_se=%.10f\n", b[["post"]], se[["post"]]))
cat(sprintf("slope_change=%.10f\nslope_change_se=%.10f\n", b[["months_since"]], se[["months_since"]]))
cat(sprintf("single_series_level_change=%.10f\nsingle_series_slope_change=%.10f\n", single[["post"]], single[["months_since"]]))
cat(sprintf("n=%d\n", nrow(d)))
