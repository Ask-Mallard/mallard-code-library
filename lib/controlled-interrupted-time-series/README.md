# Controlled interrupted time series

A single-series ITS credits everything that changes at the start date to the intervention. A control
series, exposed to the same other changes but not to the intervention, removes them.

| This fixture | Estimate | Truth |
|---|---|---|
| **level change, relative to control** | **−3.90** | −4 |
| slope change, relative to control | −0.31 | −0.3 |
| level change, intervention series alone | −5.81 | (−4, plus a shared drop of −2) |

## The model

Segmented regression of the **difference series** (intervention minus control). With both sites observed
on the same months it gives exactly the group-by-period coefficients of the full interaction model on
the stacked data (the control file checks it), and its residuals are one ordered series, so Newey-West
applies directly: lag 3, Bartlett, no pre-whitening, no adjustment.

## What it assumes

Without the intervention, the difference between the sites would have continued its pre-period trend.
Choose a control that shares the intervention site's other exposures.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
