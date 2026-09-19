# Two-by-two table: chi-square, Fisher's exact test and the odds ratio

The crude table. No model, no adjustment, nothing to misspecify.

## The default this entry exists to pin

**R and Python both apply Yates' continuity correction to a 2x2 by default, and every file here
switches it off.**

| Language | The obvious line | With the correction off |
|---|---|---|
| R | `chisq.test(tbl)` — **corrected**, 49.185 | `chisq.test(tbl, correct = FALSE)` — 49.901 |
| Python | `chi2_contingency(tab)` — **corrected**, 49.185 | `chi2_contingency(tab, correction=False)` — 49.901 |

The two conventions differ by 0.715 on this fixture, and no output says which one it used.

**Every expected count in this fixture is in the hundreds.** That is deliberate: the correction is
defended as a small-sample device, so a sparse fixture would let the disagreement be dismissed as an
edge case. It is not one. With two thousand patients the packages still differ by 0.715.

Uncorrected is what both files pin, because the correction approximates an exact conditional
test and is conservative to a fault where the expected counts are large. Where they are *not* large,
the honest answer is Fisher's exact test — computed in every file here — rather than a corrected
approximation of an asymptotic one.

## Two estimators wearing one name

`fisher.test()$estimate` in R is the **conditional maximum likelihood** odds ratio.
`fisher_exact()`'s statistic in Python is the **sample** odds ratio `ad/bc`. They are different
estimators of the same parameter and they do not agree.

So two files each written as "report the odds ratio from Fisher's test" disagree, and the
disagreement reads as a bug in one of them. Every file here computes the odds ratio **from the
table** instead, and says so where a reader will see it. `scipy.stats.contingency.odds_ratio` is the
function that matches R's, for anyone who wants the conditional estimate.

## The orientation trap

R's `table()` orders levels by value, so `0` comes before `1` and the table
arrives with the unexposed row and the no-event column first.

That leaves the **odds ratio unchanged** — swapping both rows and columns leaves `ad/bc` alone — and
silently reverses what the **risk difference** means, which then reads as the risk of *not* having
the event among the *unexposed*. A quantity that is right by luck beside one that is wrong in
silence is the worst of both. The R file states the factor levels explicitly, so a reader can see
the ordering rather than having to remember an option.

## Fisher's p is compared on the log scale

The p value here is 1.9e-12. A tolerance of 1e-4 on a number that small is satisfied by any two
answers at all, including two wrong ones, so the harness key is `neg_log10_p_fisher`.

Worth flagging: R and scipy each resolve near-ties among candidate tables with their own relative
tolerance, so this is the key most likely to expose a genuine difference between two correct
implementations rather than a mistake in either. It is included anyway. A cross-language check whose
keys are chosen to avoid the awkward ones is a check calibrated against nothing.

## What is not in this entry

**Adjustment.** This is the crude table, so it estimates a causal contrast only where the exposure
was randomised or exchangeability can be defended without covariates. Where adjustment is needed the
entry is logistic regression or the modified Poisson model, not this one with a covariate bolted on.

**Paired or matched data.** The same patients before and after, or cases matched to controls, breaks
the independence this test assumes. McNemar's test is the paired counterpart and conditional logistic
regression the matched one, which is its own entry.

**Tables bigger than 2x2.** Yates' correction does not arise there, and the trap this entry exists
for does not either.
