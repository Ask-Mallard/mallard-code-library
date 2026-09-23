# Kruskal-Wallis test with group medians

The rank-based counterpart of one-way ANOVA, for a skewed or ordinal outcome across three or more
groups.

## What it says, and what it does not

The test says whether the groups differ in distribution. It gives **no effect size**, so the group
medians are reported beside it. If specific pairs matter, estimate each pair's shift with
`mann-whitney-hodges-lehmann` and adjust for multiplicity. The test reads as a comparison of medians
only when the groups share one shape.

## Defaults

Both languages apply the tie correction and use the chi-square approximation with k − 1 degrees of
freedom, and both split the data by arm label. There is no default here that silently changes the
answer; the risk this entry guards is reporting the p-value alone.

## The fixture

Three arms of 80, log-normal with sigma 0.5 and true medians e^1.0, e^1.3 and e^1.6 (2.72, 3.67, 4.95).

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
