# Case-cohort analysis with a weighted Cox model

When measuring an exposure in a whole cohort is too expensive (a stored-sample assay, a chart
abstraction), measure it in a **random subcohort** and in **every case**. The analysis weights the
subcohort to stand for the whole cohort and uses a robust variance for the sampling.

## Why the weighting matters

| Analysis | Log hazard ratio | SE |
|---|---|---|
| **case-cohort, Lin-Ying weights** | 0.700 | 0.135 |
| case-cohort, Prentice | 0.676 | 0.137 |
| unweighted Cox on the sample (treated as the cohort) | 0.431 | 0.080 (falsely small) |
| full cohort (never available in a real study) | 0.602 | 0.080 |
| truth | 0.588 | |

## The defaults this entry pins

- `cohort.size` is the **full** cohort's size, not the sample's.
- `method = "LinYing"`; `cch`'s default is Prentice.
- `cch` returns an **unnamed** coefficient vector; the code indexes it by position and checks its
  length.
- The subcohort must be a simple random sample of the cohort for these weights. A stratified subcohort
  needs the Borgan estimators (`I.Borgan`, `II.Borgan`).

## R only

No standard Python package implements case-cohort estimators with their variance, and a hand-weighted
lifelines fit would depend on a robust variance with delayed entry that lifelines 0.28.0 gets wrong
(see `andersen-gill-recurrent-events`). Recovery is checked by fitting 100 generated fixtures in R, and
the full cohort serves as the oracle.

## The fixture

A cohort of 4000, 30% exposed, hazard 0.03 a year × 1.8 if exposed, followed up to 5 years; a 10%
random subcohort (400) plus all 643 cases: 974 rows.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | not applicable (see above) |
