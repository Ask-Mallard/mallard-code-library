# Nested case-control with risk-set sampling

When measuring exposure on a whole cohort is too expensive, take every case and, **at the moment each case
occurs**, sample a few controls from those still at risk (incidence density sampling). Conditional logistic
regression on those matched sets estimates the cohort's **hazard ratio**, with no rare-disease assumption.

| This fixture (cohort of 3,000, outcome in about 40%) | log |
|---|---|
| cohort hazard ratio, truth | 0.693 (HR 2) |
| **nested case-control, 1,230 sets of 1:4** | **0.717** (SE 0.069) |
| cumulative sampling (controls = everyone event-free at the end), mean over 100 cohorts | 0.940 |

## Keep the sets together

The analysis is conditional on the matched set (`clogit(... + strata(set))`). A person can be a control in
several sets and later a case; that is correct, not a duplicate to remove.

## Sampling from the wrong pool

Controls drawn from those who never had the event by the end of follow-up estimate an odds ratio for
cumulative risk, which departs from the hazard ratio when the outcome is common. The fixture's outcome is
common on purpose, to make that visible.

## The defaults this entry pins

R `clogit(method = "exact")`; Python `ConditionalLogit(method = "newton")`, since its default optimizer
stops short (see `case-crossover-conditional-logistic`).

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
