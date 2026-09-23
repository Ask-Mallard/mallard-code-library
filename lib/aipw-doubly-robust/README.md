# Augmented inverse-probability weighting (AIPW)

AIPW combines a propensity model and an outcome model:

    mu1 = mean( m1 + A (Y − m1) / e ),   mu0 = mean( m0 + (1 − A)(Y − m0) / (1 − e) )

It is **doubly robust**: consistent if either model is right. Weighting alone needs the propensity model;
g-computation alone needs the outcome model.

| Over 100 seeds, bias of the risk difference | AIPW | single model |
|---|---|---|
| propensity model broken (intercept only) | −0.002 | +0.125 (weighting) |
| outcome model broken (treatment only) | −0.002 | +0.125 (g-computation) |

## The standard error

From the influence function, sd(IF)/√n, which ignores that the two models were estimated. That is
valid when both are right (AIPW is then efficient); when only one is, bootstrap the whole procedure.
Packages such as `AIPW` and `tmle` add cross-fitting and flexible models around the same estimator.

## The fixture

2,500 patients; severity and a comorbidity confound treatment and a binary outcome. True marginal risk
difference 0.146; the committed estimate is 0.103 (SE 0.019), about two standard errors low, which the
calibrated tolerance allows.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
