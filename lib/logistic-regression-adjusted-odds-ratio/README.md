# Logistic regression for an adjusted odds ratio

Fits `outcome ~ exposed + covariate` by maximum likelihood with a logit link, and reads the
coefficient on `exposed` as a **log odds ratio, conditional on the covariates in the model**.
Confidence intervals are reported two ways, profile-likelihood and Wald, because the two are not
the same interval and the languages do not default to the same one.

## What this odds ratio is not

**It is not a risk ratio.** With a 23.6% event rate in this fixture the odds ratio overstates the
risk ratio substantially: standardised over the same data, the odds ratio is 1.94 and the risk
ratio is 1.66. If the question is about risk, `modified-poisson-adjusted-risk-ratio` in this
library fits the risk directly rather than converting an odds ratio afterwards.

**It is not collapsible.** This is the property that most often gets lost, so the fixture prints it
rather than describing it. From the generating model, on one dataset:

| quantity | value | what it is |
|---|---|---|
| conditional odds ratio | **2.0000** | what logistic regression estimates, and what `truth` states |
| crude odds ratio | 2.6062 | unadjusted — the gap from 2.0 is **confounding** |
| marginal odds ratio | 1.9403 | the true model standardised over the observed covariates |
| marginal risk ratio | 1.6602 | a different estimand again |

The third row is the one to sit with. It is computed from the **true** model by averaging both
potential outcomes over every row, so no confounding survives in that comparison and nothing is
estimated. It is still not 2.0. That gap is not bias and not sampling error: a conditional odds
ratio and a marginal odds ratio are different numbers, and adding a covariate to a logistic model
moves the coefficient **even when that covariate confounds nothing**. An odds ratio is therefore
not transportable between models with different covariate sets, and two correct papers reporting
different adjusted odds ratios for the same exposure are not necessarily in conflict.

This is also why Baron-Kenny mediation and product-of-coefficients decomposition are invalid for a
binary outcome, and why a marginal effect wants standardisation (g-computation) or a marginal model
rather than a conditional coefficient.

## The default this entry exists to pin

The interval, and the difference is **between languages** rather than inside one:

| call | what it returns |
|---|---|
| R `confint(fit)` | **profile likelihood** |
| R `confint.default(fit)` | **Wald** |
| Python `fit.conf_int()` | **Wald** |
| SAS `oddsratio` without `cl=pl` | **Wald** |
| Stata `logit`/`logistic` | **Wald**, with no built-in profile option |

So the obvious one-liner — the call a reader writes without thinking about it — returns a different
*kind* of interval in R than in Python, and nothing in either output says so. R is the odd one out,
and it is the one that is odd by being better: profile-likelihood intervals respect the curvature of
the likelihood and do not have to be symmetric, so they behave when the Wald interval does not.

statsmodels has no profile interval for a GLM, so `python.py` computes one: it holds the exposure
coefficient at a candidate value with an **offset**, refits everything else, and solves for where
the deviance rises by `qchisq(0.95, 1)`. That is the profile deviance itself rather than an
approximation to it.

**Measured, because asserting it would be the thing this library exists to stop.** R interpolates
the profile crossing with a spline; Python solves for the root. On identical rows the two land
within **5.9e-7** of each other on the worst bound, so the shared 1e-4 tolerance needed no
loosening and none is declared.

Both files also check that the outcome is coded 0/1 instead of trusting it. R's
`glm(family = binomial)` models the **second** level of a factor, so an outcome arriving as
`factor(c("yes","no"))` is modelled as `P(yes)` only by the accident of alphabetical order.

## The fixture

`fixture.py` writes 4000 rows with a 23.6% event rate. The outcome is drawn from a logit-link
model, so the true **conditional** odds ratio is exactly 2.0 and is a property of the data rather
than of the analysis. Exposure depends on the covariate, so the covariate genuinely confounds and
an entry about adjustment has something to adjust for.

## What the checks actually catch

Three negative controls were run against the real R output, and each fails:

| substitution | caught by | margin |
|---|---|---|
| unadjusted (crude) fit | recovery | misses truth by 0.265 against a 0.25 limit |
| probit link | both | misses by 0.327 on the exposure, 0.200 on the covariate |
| Wald bound labelled profile | **agreement only** | spread 2.3e-4 against a 1e-4 limit |

The last row is the division of labour between the two claims: swapping the interval type leaves
the point estimate untouched, so no recovery check could ever see it, and only running both
languages and comparing the bounds does.

Stated honestly: at n=4000 with 946 events the profile and Wald intervals have nearly converged,
so that 2.3x margin is thinner than the 40x the modified-Poisson entry has. The two separate when
events are few — which is exactly when the choice between them matters most, and exactly what a
well-powered fixture is least able to show.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
| SAS | not executed |
| Stata | not executed |

A working implementation of logistic regression says nothing about whether an adjusted odds ratio
answers your question, whether the covariates in the model are the right ones, or whether
conditioning on them opened a path rather than closing one.
