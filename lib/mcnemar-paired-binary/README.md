# McNemar's test for a paired binary outcome

The same people tested twice (before and after, or two tests on each person). Only the **discordant**
pairs, positive then negative or negative then positive, carry information about change.

Report the **change in the proportion positive**, (c − b) / n, with its paired Wald interval.
McNemar's chi-square and the exact binomial test on the discordant pairs are the tests; with few
discordant pairs, use the exact p-value.

## The defaults this entry pins

| Language | Default | Pinned |
|---|---|---|
| R `mcnemar.test` | continuity correction on | `correct = FALSE`; the exact test covers small counts |
| Python `statsmodels mcnemar` | `exact=True` | the chi-square as `exact=False, correction=False`; the exact p-value named |

A two-sample chi-square on the before and after columns treats them as independent groups. On this
fixture it gives 4.68 against McNemar's 13.13.

## The fixture

500 people; 30% positive before; after, 80% stay positive and 20% of negatives turn positive. True
change +0.08. Observed: 23 positive→negative, 55 negative→positive, change 0.064.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
