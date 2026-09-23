# Overlap weights for the effect in the overlap population (ATO)

Each treated patient is weighted by 1 − e and each untreated patient by e, where e is the propensity
score (Li, Morgan and Zaslavsky 2018). Patients whose treatment was close to a coin toss count most.
The weights are bounded by 1, so there is nothing to trim, and with a logistic propensity model they
balance every covariate's mean **exactly**.

| This fixture | Value |
|---|---|
| ATE (everyone), truth | 2.00 |
| **ATO (overlap population), truth** | **2.25** |
| estimate | 2.47 (SE 0.13) |
| weighted SMD of severity | 0 |

Report it as the **ATO**: when the effect varies, it is a different population from the ATE.

## The standard error

The propensity model and the two weighted means are solved as one stack of estimating equations, and
the sandwich covers the whole stack, so the estimation of the score is included (as R's `PSweight`
does). The robust SE that treats the weights as known is reported beside it (0.186 here); for the ATO it
is not guaranteed to be conservative. Both languages write the stack out; the control file rebuilds it
with numerical derivatives.

## The fixture

3,000 patients; severity and a comorbidity raise treatment and outcome; the effect is 2 + 2 × severity.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
