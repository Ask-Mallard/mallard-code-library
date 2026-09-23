# Adjusted risk difference from an identity-link binomial model

A risk difference adjusted for covariates. With an **identity link**, the treatment coefficient of a
binomial model is a risk difference on the probability scale. The default logit link would give a log
odds ratio instead (0.83 on this fixture, against a risk difference of 0.095).

## The defaults this entry pins

| Language | Default | Pinned |
|---|---|---|
| R `glm` | logit link; start values from the logit scale | `binomial(link = "identity")`, start values, `epsilon = 1e-12` |
| Python `GLM` | logit link | `Binomial(link=Identity())`, start values, `tol=1e-12` |

At the default convergence settings the two languages stopped **1e-5 apart** on this fixture; with the
tolerance pinned they agree to 2e-7.

## When it fails, and what to do

Identity-link models can fail to converge, or produce risks outside [0, 1], when risks approach 0 or 1.
Both files check that fitted risks stay inside (0, 1). If a model will not converge, estimate the risk
difference by **marginal standardization (g-computation)** from a logistic model, or a linear
probability model with robust errors. Dropping covariates to force convergence changes the estimand.

## The fixture

1000 people; risk = 0.10 + 0.08 × treated + 0.004 × (age − 60), ages 40 to 80, so risks run from
0.02 to 0.26. True risk difference 0.08. Observed 0.095.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
