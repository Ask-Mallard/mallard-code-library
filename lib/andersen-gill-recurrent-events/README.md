# Andersen-Gill model for recurrent events

For events that can happen more than once: admissions, falls, infections, exacerbations. The
Andersen-Gill model uses **every** event, not just the first, and gives a population-average rate
ratio. Because some people have many events and others none, the events of one person are correlated,
and the variance must be **clustered on the person**.

| Standard error of the log rate ratio | This fixture |
|---|---|
| model-based (naive): treats every admission as independent | 0.0557 (too small) |
| **robust (Lin-Wei), clustered on the person** | **0.0701** |
| lifelines 0.28.0, `robust=True` with `entry_col` | 0.0472 (wrong: smaller than naive) |

## A robust variance that is not

lifelines' `CoxPHFitter(entry_col=..., cluster_col=..., robust=True)` reported a standard error smaller
than the naive one on this fixture. A clustered sandwich cannot be smaller here, where people differ
widely in how often they are admitted. A from-scratch Lin-Wei sandwich (each row's score residual,
summed by person) and R's `coxph(..., cluster = id)` both give 0.0701. The Python file takes the point
estimate from lifelines and writes the sandwich out; a control keeps checking the lifelines value.

## Alternatives

Time to the first event discards most events. A negative binomial on the total count
(`negative-binomial-rate-ratio`) is often simpler when event timing does not matter. A shared frailty
model estimates a conditional (person-specific) ratio instead of this population-average one.

## Learn page

The Learn page this entry names, `analysis.recurrent_events`, is planned but not yet published.

## The fixture

1000 people followed 1 to 3 years; admissions at 0.8 a year × a gamma frailty (variance 0.5) × 0.7 if
treated. 1335 admissions. True log rate ratio −0.357; observed −0.355.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
