# Zero-inflated Poisson regression for counts with structural zeros

Some people can never have the event counted: they use another hospital, or are not at risk. Those are
**structural zeros**, distinct from people at risk who happened to have none. The model has two parts:
a logistic model for being a structural zero, and a Poisson model for the count among those at risk.

## What the rate ratio means

The count part's rate ratio applies **among those at risk**, not to the whole population. Say so in the
report.

## When to use it

A zero-inflated model needs a mechanism that produces structural zeros. "There are many zeros" is not
enough: a negative binomial often fits excess zeros just as well (`negative-binomial-rate-ratio`).

## The defaults this entry pins

| Language | Default | Pinned |
|---|---|---|
| R `zeroinfl` | without `\|`, the SAME covariates in both parts | `visits ~ treated \| 1`: intercept-only zero part |
| Python `ZeroInflatedPoisson` | intercept-only zero part (the opposite of pscl) | given explicitly |
| Python `ZeroInflatedPoisson` | `bse` from an analytic Hessian | standard errors from the numerical Hessian of the log-likelihood |

## A standard-error finding

statsmodels 0.14.2's `ZeroInflatedPoisson` reported a standard error of **0.0752** for the rate ratio
here; R's pscl reported **0.0761**. An independent numerical Hessian of a written-out log-likelihood
(whose value matched statsmodels' to 1e-13) gave 0.0761. statsmodels' `bse` understated the uncertainty,
by 10% on the inflation term (0.0956 against 0.1060). The Python file uses the log-likelihood's
curvature, and a control keeps checking the discrepancy.

## The fixture

800 people; 30% structural zeros; otherwise Poisson with mean 2 (control) or 2 × e^−0.5 (treated). True
log rate ratio −0.5. Observed −0.44.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
