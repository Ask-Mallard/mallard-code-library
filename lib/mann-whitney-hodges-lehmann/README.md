# Mann-Whitney test with the Hodges-Lehmann shift

A rank comparison of two independent groups when the outcome is skewed or ordinal. Report the
**Hodges-Lehmann shift** with its interval, not the p-value alone: the p-value says the groups
differ, and the shift says by how much.

## What the shift is, and is not

The Hodges-Lehmann estimate is the **median of all treated-minus-control pairwise differences**. It is
not the difference between the two medians (on this fixture 2.49 against a shift of 2.43), and the
Mann-Whitney test compares medians only when the two groups have the same shape. When the shapes
differ, the test compares the probability that a random treated value exceeds a random control value.

## The defaults this entry pins

| Language | Default | Consequence |
|---|---|---|
| R `wilcox.test` | `conf.int = FALSE` | no estimate and no interval, only a p-value |
| R `wilcox.test` | exact only below 50 per group with no ties | a normal approximation with continuity correction otherwise |
| Python `mannwhitneyu` | `method="auto"` | a normal approximation above 8 per group |
| Python | no Hodges-Lehmann interval | built here from the exact null distribution of U |

Treated is passed first in both languages; the sign of the shift follows argument order.

## Ties

This fixture has none. With ties (common for ordinal outcomes and rounded data), the exact
distribution does not apply. R warns and falls back to a normal approximation; scipy's
`method="exact"` does not account for ties at all, so use `method="asymptotic"` there. The
exact-order-statistic interval built in `python.py` assumes no ties too. Report which was used.

## The fixture

40 per group, one log-normal shape, the treated group shifted up by exactly 2. Observed shift 2.43.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
