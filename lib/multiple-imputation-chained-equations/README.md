# Multiple imputation by chained equations, pooled by Rubin's rules

A confounder (baseline severity) is missing for a third of patients, more often when the outcome is high.
That is **missing at random** given the outcome. Dropping the incomplete rows biases the adjusted effect;
imputing severity 20 times, with the outcome in the imputation model, and pooling by Rubin's rules does
not.

| This fixture | Treatment effect (truth 1.0) |
|---|---|
| **multiple imputation** | **1.00** (df 214, fraction of missing information 0.26) |
| complete-case analysis | 0.82 |

Over 100 simulated datasets MI averages 0.99 and complete-case analysis 0.74.

## Name the assumption first

MI is valid under MAR. Under missing not at random it is not enough on its own: add a delta-adjusted
sensitivity analysis (`delta-adjusted-tipping-point`). MI is not better in itself; it is right when its
assumption is.

## The defaults this entry pins

`mice` with **m = 20** (its default is 5), maxit = 10, method `norm` (Bayesian linear regression) for the
incomplete covariate, **the outcome in the imputation model**, a fixed seed, and `pool()` (Rubin's rules
with Barnard-Rubin degrees of freedom). Report the fraction of missing information.

## R only

Imputation draws are random and R's and Python's generators and MICE implementations differ, so the two
cannot be compared draw for draw. The control file recomputes Rubin's rules and the Barnard-Rubin df from
mice's pooled components.

## Verification

| Engine | Status |
|---|---|
| R (mice 3.19.0) | executed in CI |
| Python | not applicable (see above) |
