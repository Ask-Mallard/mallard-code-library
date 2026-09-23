# Penalized prediction model development with bootstrap optimism correction

A model fitted and assessed on the same patients looks better than it will in new ones. Two remedies,
used together:

1. **Shrink at development.** Ridge (L2-penalized) logistic regression, the penalty chosen by 5-fold
   cross-validation.
2. **Correct the apparent performance by the bootstrap, repeating every step** (Harrell): in each
   resample, choose the penalty again, fit, and record the AUC and calibration slope in the resample and
   in the original data. The average difference is the optimism.

| This fixture (300 patients, 68 events, 10 candidate predictors) | AUC | Calibration slope |
|---|---|---|
| apparent | 0.753 | 1.23 |
| optimism-corrected | 0.713 | 0.99 |
| external cohort (5,000 patients) | 0.746 | 1.22 |

Over 100 simulated development samples the correction is **unbiased**: apparent AUC exceeds the
external AUC by 0.033 on average, the corrected one by −0.0015 (slope: +0.20 against +0.0006). In any
single sample it is **imprecise**, as this fixture shows. Sixty-eight events for 10 candidate predictors
is below the Riley et al. (2019) criteria; the fixture is small on purpose, to make the optimism visible.

## The defaults this entry pins

The ridge fit is written out in both languages (Newton-Raphson; objective loglik − (λ/2)Σβ² on
standardized predictors, intercept unpenalized), because glmnet and scikit-learn scale the penalty
differently. The folds and the 50 bootstrap resamples come from the fixture (`bootstrap.csv`) so both
languages use the same ones; use 200 or more in practice.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
