# Fisher's exact test with the conditional odds ratio

An exact test for a 2x2 table, valid at any size and the usual choice when expected cell counts are
small. Report the **odds ratio with its interval**, not the p-value alone.

## Which odds ratio

| Quantity | What reports it |
|---|---|
| **Conditional MLE** (conditioning on both margins) | R `fisher.test`; scipy `contingency.odds_ratio(kind="conditional")` |
| Sample odds ratio ad/bc | scipy `fisher_exact`'s statistic |

They differ, most in small tables (5.62 against 5.67 here). Say which one a report gives.

## The defaults this entry pins

| Language | Trap | Pinned |
|---|---|---|
| R | `table(exposed, outcome)` sorts both to 0 first; with one of the two flipped the odds ratio inverts | the 2x2 built explicitly, exposed row first, outcome column first |
| Python | `fisher_exact`'s statistic is the sample odds ratio | the odds ratio from `odds_ratio(kind="conditional")`; `fisher_exact` supplies only the p-value |

## Agreement

The estimate and exact bounds are roots of equations that R and scipy solve with different
root-finders, so they agree closely rather than exactly (the upper bound to about 1e-4 on the log
scale). `expected.json` records the measured gaps and the tolerance.

## The fixture

120 per group; risks 0.5 and 0.2, a true odds ratio of 4. At 30 and 60 per group the recovery check
could not tell an odds ratio of 4 from 1, so the fixture is larger than a "small table" example.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
