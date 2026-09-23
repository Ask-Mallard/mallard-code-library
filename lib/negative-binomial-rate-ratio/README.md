# Negative binomial regression for overdispersed counts

Event counts (falls, admissions, exacerbations) usually vary more than a Poisson model allows, because
some people are prone to events and some are not. A Poisson model then gets the rate ratio right and
its standard error **too small**. The negative binomial (NB2) model estimates the extra variance
through a dispersion parameter, theta.

| Model | SE of the log rate ratio on this fixture |
|---|---|
| Poisson | 0.061 (too small) |
| NegativeBinomial family at its default alpha = 1 | 0.090 (an assumed dispersion) |
| **Negative binomial, dispersion estimated** | **0.079** |

## The defaults this entry pins

| Language | Trap | Pinned |
|---|---|---|
| both | no offset: the model compares counts, and longer follow-up looks worse | `offset(log(follow-up))` |
| Python | `sm.families.NegativeBinomial()` fixes alpha = 1 | dispersion estimated by the discrete `NegativeBinomial`, then fixed in a GLM, which is what `glm.nb` does |
| both | default convergence stopped theta 7.5e-5 apart | tolerances tightened to 1e-12 |

## The fixture

1000 people followed 0.5 to 3 years; person-level rates are gamma-distributed around 0.8 a year
(control) or 0.56 (treated), shape 2. True rate ratio 0.7, theta 2. Observed log rate ratio −0.362,
theta 1.75.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
