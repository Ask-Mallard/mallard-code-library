# Paired t-test for a within-person change

The same people measured twice, before and after. Each person is their own control, so the analysis
is a one-sample t on the **differences**, and the result is the mean change with its interval.

## Why pairing matters

People differ from each other far more than each person changes. A two-sample (Welch) test on the
same two columns counts that between-person spread as noise: on this fixture its interval is
**4.2 times wider** for the same data.

## The defaults this entry pins

| Language | Trap | Pinned |
|---|---|---|
| R `t.test` | argument order sets the sign of the change; the formula interface no longer accepts `paired = TRUE` (R 4.4) | `t.test(after, before, paired = TRUE)` |
| Python | `ttest_ind` is the unpaired test | `ttest_rel(after, before)` |

## What it assumes, and what it cannot say

The **differences** should be roughly normal; the raw values need not be. With skewed or outlying
differences, use `wilcoxon-signed-rank-hodges-lehmann`.

A change in one group is not the effect of whatever happened between the measurements. Regression to
the mean, secular trends and repeated testing all produce change without an intervention. A causal
claim needs a comparison group.

## The fixture

60 people; before around 150 (SD 15), after = before − 5 + noise (SD 6). True mean change −5,
true SD of the change 6. Observed: −5.14.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
