# Multinomial logistic regression for an unordered categorical outcome

For an outcome with categories that have **no order**: discharge destination, tumour subtype, reason
for referral. One log relative-risk ratio per non-reference category, each comparing that category
with the reference.

## The reference category decides what every coefficient means

With home as the reference, the treatment coefficient for nursing is −0.89 here: treatment shifts
patients from nursing toward home. With rehab as the reference, the same data give −1.47 for nursing,
because that coefficient now compares nursing with rehab. Choose the reference on purpose and name it.

## The defaults this entry pins

| Language | Default | Pinned |
|---|---|---|
| R `factor()` | levels in alphabetical order | levels written out, home first |
| R `multinom` | 100 iterations, loose tolerance: stopped 2.45e-4 short here | `maxit = 1000, reltol = 1e-14` |
| Python `MNLogit` | lowest integer code is the reference | codes assigned explicitly, home = 0 |

## Nominal or ordinal?

If the categories do have an order, a proportional-odds model uses it and is more efficient
(`ordinal-proportional-odds`).

## Learn page

The Learn page this entry names, `analysis.multinomial_regression`, is planned but not yet published.

## The fixture

2000 patients; log(P(rehab)/P(home)) = −0.5 + 0.6 × treated, log(P(nursing)/P(home)) = −1.0 − 0.8 ×
treated. Observed 0.58 and −0.89.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
