# An incidence rate with an exact Poisson interval

Events per person-year of follow-up in a cohort, with a 95% interval. When people are followed for
different lengths of time there is no single "risk"; the rate, with **person-time** as its
denominator, is what the data estimate.

## The default this entry pins

| Language | Default | Consequence of relying on it |
|---|---|---|
| R `poisson.test` | `T = 1` | a valid interval for the event **count**, not the rate, with no warning |
| Python `confint_poisson` | exposure of 1; about a dozen methods | the count's interval, or a method you did not mean |

The interval is the central exact (Garwood) interval: `method="exact-c"` in statsmodels, the default
in R's `poisson.test`.

## When not to use it

The exact interval assumes events are **Poisson**: independent, at one rate for everyone. Cohorts
where some people have repeated events (falls, admissions, infections) are usually overdispersed,
and this interval is then too narrow. Use a negative binomial model or a rate regression with a
robust (sandwich) interval instead (`poisson-rate-regression-overdispersed`).

## The fixture

600 people followed for 0.5 to 5 years each, events arriving as a Poisson process at 0.08 per
person-year. Follow-up varies on purpose: with equal follow-up, dividing by people instead of
person-time would give the same answer after rescaling. Observed: 146 events over 1,617 person-years,
0.0903 per person-year.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |

`harness/descriptive_examples.test.py` checks the bounds against Garwood's chi-square formula and
measured the recovery tolerance over 100 seeds.
