# One-way ANOVA with Tukey pairwise intervals

Comparing a continuous outcome across three or more groups. The **F test** asks whether any group
mean differs; it does not say which. The **pairwise differences with Tukey-adjusted intervals** are
the result a reader needs: with three comparisons, the family of three intervals keeps 95% coverage
together.

## The default this entry pins

The arms are coded 1, 2 and 3. A model that reads that column as a **number** fits one straight line
(1 degree of freedom) instead of three group means (2 degrees of freedom), and prints an ordinary
ANOVA table. On this fixture the F statistic goes from 94.1 to 176.6.

| Language | Pinned |
|---|---|
| R | `aov(y ~ factor(arm))`, then `TukeyHSD` |
| Python | groups split by label, then `f_oneway` and `tukey_hsd` |

## What it assumes

Roughly normal outcomes with similar SDs in each group. With clearly unequal SDs prefer Welch's
ANOVA (`oneway.test`) and Games-Howell comparisons; with skewed outcomes, `kruskal-wallis`.

## The fixture

Three arms of 100, true means 10, 12 and 15, common SD 2.5. The means are unequally spaced (2 then 3)
so a straight-line fit is visibly wrong.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |

The Tukey bounds agree to 3e-10: R and scipy integrate the studentized range distribution
independently.
