# Bayesian logistic regression with stated priors

The same model as a frequentist logistic regression, with priors on the coefficients and the answer
reported as a posterior distribution.

| This fixture (1,000 patients) | Value |
|---|---|
| true treatment log OR | 0.60 |
| **posterior mean (SD)** | **0.585 (0.137)** |
| 95% credible interval | 0.317 to 0.858 |
| P(log OR > 0) | 1.000 |
| R-hat, effective sample size | 1.000, about 2,000 |
| maximum-likelihood estimate | 0.59 |

## What to report

- **The priors, and why.** Here normal(0, 2.5) on each coefficient and normal(0, 5) on the
  intercept: weakly informative on the log-odds scale, ruling out implausibly huge odds ratios without
  favouring a direction. Add a sensitivity analysis
  under a sceptical prior when the conclusion could depend on the prior.
- **The posterior**: mean and SD, a credible interval, and the probability of benefit (or of a clinically
  important effect), not a p-value.
- **Sampler diagnostics**: several chains, R-hat near 1, an effective sample size in the thousands.

## The sampler

Random-walk Metropolis written out in R (4 chains, 5,000 draws each after warm-up) so every step is
visible and nothing needs compiling. In practice use Stan (`rstanarm::stan_glm` or `brms`), which samples
far more efficiently; the model and priors are the same.

## R only, truth check only

The Bayesian entry is checked against its truth only: Markov chains cannot be
matched across languages. The control file compares the posterior with a Laplace approximation computed
in Python, and over 100 simulated trials the 95% credible interval covered the true effect 95% of the
time.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | not applicable (see above) |
