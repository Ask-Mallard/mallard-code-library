# Table 1: baseline characteristics with standardized differences

The table every study opens with: who was included, overall and by group.

## Summarize each variable by its type, not its position

| Variable type | Report | Example here |
|---|---|---|
| continuous, roughly symmetric | mean (SD) | age |
| continuous, skewed | median [Q1, Q3] | length of stay |
| binary | n (%) | diabetes |
| categorical | n (%) for each level | smoking: never, former, current |

## Compare groups with standardized differences, not p-values

A p-value in Table 1 tests a difference the reader can already see, and it depends on sample size.
In a randomized trial it tests a null hypothesis that is true by design, and CONSORT advises against
it. The **standardized mean difference** (SMD) describes the size of an imbalance on a common scale.
An absolute SMD above about 0.1 is the usual flag for imbalance worth adjusting for.

The SMD pools the group SDs as sqrt((s1² + s0²) / 2), and a binary variable uses p(1 − p) as its
variance (Austin 2009).

## The defaults this entry pins

| Language | Default | Consequence of relying on it |
|---|---|---|
| numpy `np.std` | `ddof=0`, divides by n | a smaller SD than R's `sd` and pandas' `std`, which divide by n − 1 |
| quantiles | nine definitions in use; SAS and SPSS default to others | quartiles that differ between packages on the same data |

Both are pinned: `ddof=1`, and R's `type = 7` with pandas' `interpolation="linear"`.

## Missing values

This fixture is complete. With missing values, report the number missing for each variable and say
whether percentages use the non-missing count.

## The fixture

600 people, exposure assigned with probability 0.5. Age differs by group (60 against 65, SD 11), so
the true age SMD is 0.45. Diabetes is 25% against 35%, a true SMD of 0.22. Length of stay is
log-normal with a median of 3.32 days and the same in both groups.

## Learn page

The Learn page this entry names, `analysis.descriptive_statistics`, is planned but not yet published.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |

`harness/descriptive_examples.test.py` checks every value against Python's `statistics` module, shows
each pinned default changes a number when dropped, and measured every recovery tolerance over
100 seeds.
