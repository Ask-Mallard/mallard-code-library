# Restricted cubic spline for a nonlinear exposure

A smooth curve for a continuous exposure such as age, instead of a straight line (usually wrong) or
categories (which throw information away and put artificial steps at the cut points). The curve is
cubic between the outer knots and **linear beyond them**, which keeps the tails from swinging wildly.

## Report the curve, not the coefficients

Spline coefficients depend on the basis and mean nothing on their own. Report predicted values at chosen
exposures, and contrasts between them with intervals: "systolic pressure is 22.6 mmHg higher at 70 than
at 50 (95% ...)". Here R and Python use **different bases** (B-spline and truncated-power) and agree to
ten digits on every prediction, because both span the same curves.

## The default this entry pins

`ns(age, df = 4)` places knots at the data's quantiles. Software computes quantiles differently, so the
same call can fit different curves in different packages: quantile-placed knots moved the 70 vs 50
contrast from 22.64 to 21.53 on this fixture. Both files write the knots out (30, 45, 60, 70, 80).
Harrell's usual choice is knots at fixed quantiles of the exposure; choose them from the design, then
write them down.

## The fixture

600 adults aged 20 to 90; the true mean pressure is itself a restricted cubic spline with the same
knots, rising from about 115 to 150 mmHg around age 58, so the true curve has no approximation error.
True 70 vs 50 contrast 21.23; observed 22.64.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
