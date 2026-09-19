# Diagnostic accuracy with Wilson intervals

Sensitivity and specificity against a reference standard. Each is a **single proportion**, which is
precisely the case where a Wilson or Clopper-Pearson interval is the right choice: the
normal-approximation (Wald) interval degrades near 0 and 1 and can run outside [0, 1], which for an
accuracy measure is visibly wrong on the page.

**That licence does not extend anywhere else.** Wilson and Clopper-Pearson belong to a single
proportion or rate. An odds ratio, a hazard ratio or any adjusted model estimate takes a
profile-likelihood or Wald interval, and a bootstrapped quantity takes a percentile or
bias-corrected one.

## Why the intervals are the interesting part

Both languages default to something other than Wilson:

| Language | Default | Consequence of omitting the option |
|---|---|---|
| R `binom.confint` | `methods = "all"` | returns **eleven** intervals; you report whichever row you indexed |
| Python `proportion_confint` | `"normal"` | the Wald interval, which can leave [0, 1] |

R and Python use **independent implementations** of the Wilson formula rather than the same
arithmetic typed twice, so their bounds agreeing is evidence rather than a tautology.

## The fixture

2000 rows. Disease status is drawn first, then the test result conditional on it, so sensitivity
and specificity are properties of the generator.

**Prevalence is 0.30, not 0.50, deliberately.** Sensitivity and specificity do not depend on
prevalence, so at 0.50 an implementation quietly computing *predictive values* instead would still
land near the truth and pass. At 0.30 the two separate clearly.

Observed in the generated data: sensitivity **0.8507** against a true 0.85, specificity **0.9012**
against a true 0.90, prevalence 0.3015.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
