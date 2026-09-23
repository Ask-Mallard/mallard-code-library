# Interrupted time series: segmented regression with Newey-West errors

A single monthly series, before and after an intervention. Segmented regression estimates the
pre-intervention trend, a **level change** at the start, and a **slope change** after it. Report both.

| This fixture | Estimate | Truth |
|---|---|---|
| level change | −5.22 | −4 |
| slope change (a month) | −0.33 | −0.3 |

## Autocorrelation

Neighbouring months are correlated (lag-1 residual autocorrelation 0.56 here; OLS SE of the level change 0.33 against Newey-West 0.57), and ordinary least squares
standard errors are then too small. The Newey-West (HAC) covariance allows correlation up to a stated
lag. **Pinned**: lag 3, Bartlett weights, no pre-whitening, no small-sample adjustment.
`sandwich::NeweyWest`'s own defaults pre-whiten and choose the lag from the data, so name the options.
Prais-Winsten or ARIMA models are alternatives that model the correlation instead.

## What this fixture does not have

A seasonal cycle. A real monthly series usually needs seasonal terms. A single series also cannot rule
out another change at the same date; see `controlled-interrupted-time-series`.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |

The control file rebuilds the Newey-West covariance by hand.
