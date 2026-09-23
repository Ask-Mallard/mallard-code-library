# Cochran-Armitage test for trend in a proportion

Whether a binary outcome becomes steadily more (or less) common across ordered groups: dose levels,
age bands, exposure categories. It spends 1 degree of freedom on the trend, where a chi-square on the
whole table spends k − 1, so it is more powerful when the trend is real and blind to a pattern that
rises and then falls. Report the **slope** (change in proportion per score step) beside the test.

## The defaults this entry pins

| Language | Default | Pinned |
|---|---|---|
| R `prop.trend.test` | `score = 1, 2, ..., k` | `score =` the actual dose values |
| Python | no package reports this statistic | written out, dividing by N; the linear-by-linear test divides by N − 1 |

The test is only as meaningful as its scores. With unequally spaced doses, equal default scores test a
different trend.

## The fixture

Four dose groups (scores 0 to 3) of 150, proportions 0.10, 0.16, 0.22, 0.28: exactly linear, slope
0.06. Observed slope 0.051.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
