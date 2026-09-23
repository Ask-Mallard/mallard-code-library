# Case-crossover design by conditional logistic regression

Only cases are sampled. Each patient's exposure just before the event (the hazard window) is compared
with their own exposure at other times (referent windows), so everything that does not change over those
windows is removed by design. It suits transient exposures with acute effects.

| This fixture | Value |
|---|---|
| truth | odds ratio 2.50 (log 0.916) |
| **conditional logistic** | **log OR 0.899** (SE 0.111) |
| ordinary logistic, matching ignored | log OR 0.744 |
| informative patients (exposure varied across windows) | 385 of 600 |

## The defaults this entry pins

Conditional logistic regression stratified on the patient; three **time-stratified** referent windows
(fixed strata avoid the overlap bias of windows chosen relative to the event date). R: `clogit(method =
"exact")`. Python: `ConditionalLogit` with **`method="newton"`**: its default optimizer stopped 5.4e-5
short of the maximum here.

## What still biases it

A trend in exposure within the referent stratum, and exposures whose effect is not transient.

## Verification

| Engine | Status |
|---|---|
| R (survival::clogit) | executed in CI |
| Python (statsmodels ConditionalLogit) | executed in CI |

The control file finds the exact conditional maximum by Newton's method.
