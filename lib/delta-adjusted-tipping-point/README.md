# Delta-adjusted tipping-point analysis for missing outcomes

A trial's primary analysis of incomplete outcomes usually assumes **missing at random**, which cannot be
checked. A tipping-point analysis asks how far from MAR the missing outcomes would have to be to change
the conclusion.

## How it works

1. Within each arm, model the outcome on baseline covariates among patients with an outcome, and average
   observed and predicted outcomes (the MAR estimate).
2. Add a shift **delta** to the predicted outcomes of the **treated** patients whose outcome is missing
   (delta < 0: they did worse than MAR predicts).
3. The **tipping point** is the delta at which the lower 95% limit of the effect reaches zero.

| This fixture | Value |
|---|---|
| effect under MAR (truth 1.5) | 1.29 (95% CI 0.88 to 1.69) |
| complete-case difference | 0.73 |
| **tipping point** | **delta = −2.02** (treated dropouts 2.0 units worse than predicted) |

Clinicians then judge whether a shift that large is plausible. Shifting both arms, or a grid of
(treated, control) shifts, are common variants.

## The standard error

From the stacked estimating equations (each arm's regression and mean), written out in both languages,
so the estimate is deterministic. Multiple imputation with a delta added to the imputed values gives the
same idea with Monte Carlo error.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
