# G-computation (marginal standardization)

Fit an outcome model, predict every patient's risk as if treated and as if untreated, and average. The
result is a **population** risk difference or risk ratio, adjusted for the confounders in the model.

## The coefficient is not the answer

`exp(coef)` from a logistic model is a **conditional** odds ratio. The marginal odds ratio is smaller
even with no confounding at all (non-collapsibility):

| This fixture | Truth | Estimate |
|---|---|---|
| **risk difference** | **0.146** | 0.154 (SE 0.022) |
| log risk ratio | 0.363 | 0.375 |
| log marginal odds ratio | 0.609 | 0.641 |
| log conditional odds ratio | 0.700 | 0.730 |
| crude risk difference | | 0.274 |

## The standard error

The model and the two standardized risks are one stack of estimating equations, and the sandwich
covers all of it (as R's `stdReg`), including the sampling of the covariates. The delta method on the
model's covariance alone (the default of `marginaleffects::avg_comparisons`) treats the covariates as
fixed; it is reported beside it and is slightly smaller here.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |

The control file rebuilds the sandwich with numerical derivatives.
