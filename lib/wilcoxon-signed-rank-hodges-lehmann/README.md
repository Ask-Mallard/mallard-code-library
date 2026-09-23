# Wilcoxon signed-rank test with the Hodges-Lehmann pseudo-median change

For paired data (the same people before and after) whose changes are **symmetric but not normal**:
heavy tails, a few large changes. Report the **pseudo-median change** with its interval: the median
of all pairwise averages of the changes (Walsh averages), which under symmetry estimates the centre
of the change distribution.

## The defaults this entry pins

| Language | Default | Consequence |
|---|---|---|
| R `wilcox.test` | `conf.int = FALSE` | no estimate, no interval |
| R `wilcox.test` | exact only without zeros or ties and below 50 pairs | a normal approximation otherwise |
| Python `scipy.stats.wilcoxon` | two-sided statistic is **min(W+, W−)** | not R's V (W+). On this fixture 228 against R's larger V |
| Python | no pseudo-median or interval | built from the exact null distribution of the signed-rank statistic |

## What it assumes

Changes symmetric about their centre. With clearly skewed changes, neither the test nor the
pseudo-median answers a question about a typical change. **Zero changes** need a stated rule (R
drops them; Pratt's method keeps them); this fixture has none.

## The fixture

45 people; changes drawn from a Laplace distribution centred on +3 (symmetric, heavy-tailed). The
changes are mostly positive on purpose: only then do scipy's statistic and R's V differ, which is
what makes the pinned convention checkable. Observed pseudo-median 2.55.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
