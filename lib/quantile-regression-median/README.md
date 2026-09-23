# Median regression for a skewed outcome

For a skewed outcome such as length of stay, the question is often about a **typical** patient: "the
median stay is 1.5 days longer". Median (quantile) regression estimates that with covariate adjustment,
which a rank test cannot.

## Median effect is not mean effect

When the spread differs between groups, the two answer different questions. On this fixture the treated
vary more, and:

| Method | Estimates | Truth | This fixture |
|---|---|---|---|
| least squares | the mean effect | 1.80 | 1.66 |
| **median regression** | the median effect | 1.50 | 1.27 |

## The defaults this entry pins

| Language | Default | Pinned |
|---|---|---|
| R `summary.rq` | `se = "rank"` below 1001 rows: an interval, no standard error | `se = "nid"`, the Hendricks-Koenker local sandwich, Hall-Sheather bandwidth |
| R `summary.rq` | `"iid"` assumes one error distribution for all | not used: the spread differs by arm here |
| Python `QuantReg` | iterative approximation; kernel-sandwich standard errors unlike R's | the exact linear program (HiGHS) and the Hendricks-Koenker sandwich written out |

## An odd number of rows

With an even count the sample median is an interval, and the regression can have several optimal
solutions (R warned "Solution may be nonunique" at 400 rows). The fixture has 401 rows, so both exact
solvers return the same unique solution and agree to ten digits.

## Learn page

The Learn page this entry names, `analysis.quantile_regression`, is planned but not yet published.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
