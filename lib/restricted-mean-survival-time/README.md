# Restricted mean survival time (RMST)

The area under each Kaplan-Meier curve up to a horizon tau: the mean event-free time over those tau
years. The difference between arms ("treated patients lived 0.44 years longer on average over 4 years")
is in units a patient understands and needs **no proportional-hazards assumption**.

## The horizon

Choose tau before seeing the data, within follow-up in both arms. R's `print.survfit` reports no
restricted mean by default, and `rmean = "common"` uses the largest observed time, which moves with the
data. Both files state tau = 4 years.

## A variance that is not a standard error

lifelines' `restricted_mean_survival_time(..., return_variance=True)` returns the variance of **each
patient's** restricted survival time, min(T, tau): the spread of the data. It is not the sampling
variance of the RMST estimate.

| Quantity | Treated arm, this fixture |
|---|---|
| lifelines `return_variance`, square-rooted | **1.24 years** |
| standard error of the RMST estimate | **0.047 years** |

Used as a standard error, the first gives an interval from −3.2 to 4.1 years for a difference of 0.44.
The Python file computes the estimate's own (Greenwood-type) variance, which matches R's `se(rmean)` to
ten digits.

## The fixture

700 per arm; exponential hazards 0.20 (control) and 0.12 (treated); follow-up ends between 3 and 6
years. True RMST to 4 years: 2.753 and 3.177, a difference of 0.423. Observed 0.437.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
