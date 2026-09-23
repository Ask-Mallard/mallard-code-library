# Mixed-effects Poisson regression for clustered counts

Infection counts per patient, with length of stay as the exposure, and patients clustered in wards. A
random intercept per ward carries the between-ward variation; a Poisson model that ignores it gives a
standard error that is too small (0.096 against 0.127 here).

## Conditional or population-average

With a log link and a normal random intercept, the ward-conditional rate ratio **equals** the
population-average one: the rate ratio is collapsible. A logistic GLMM's odds ratio is not
(`mixed-effects-logistic-clustered`).

## The defaults this entry pins

`offset(log(days))` (a rate per patient-day, not a count per patient), `family = poisson`, and the
Laplace approximation (`nAGQ = 1`), lme4's default, named. Adaptive quadrature with 10 points agrees to
2e-5 on this fixture.

## R only

Python's standard stack has no frequentist Poisson GLMM; statsmodels' `PoissonBayesMixedGLM` is a
variational Bayes approximation, a different estimator. Recovery is checked by fitting 100 generated
fixtures in R.

## The fixture

40 wards (20 per arm) of 40 patients, 2 to 20 days each; rate 0.03 a patient-day × a ward multiplier
exp(N(0, 0.25)) × 0.6 in the intervention arm. True log rate ratio −0.511; observed −0.467.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | not applicable (see above) |
