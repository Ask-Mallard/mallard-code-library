# Propensity-score matching for the effect in the treated (ATT)

Each treated patient is paired with the untreated patient whose propensity score is closest. Matching
estimates the effect **in the treated**, not in everyone: when the effect varies with a confounder and
the treated differ from the untreated, those are different numbers.

| This fixture | Value |
|---|---|
| ATE (everyone), truth | 2.00 |
| **ATT (the treated), truth** | **2.68** |
| matched estimate | 2.58 (SE 0.12) |
| naive difference, unmatched | 3.43 |
| SMD of severity, before → after | 0.50 → 0.06 |

## The defaults this entry pins

MatchIt's defaults are 1:1 nearest neighbour on the score with **no caliper**. This entry names every
choice: the **logit** of the score as the distance (`link = "linear.logit"`), a **caliper of 0.2 SD** of
that logit (Austin 2011), treated patients matched **from the largest score down**, without
replacement, estimand **ATT**. A different order or caliper gives a different matched set.

The standard error clusters on the matched pair (`sandwich::vcovCL(cluster = ~subclass)`), MatchIt's
recommendation. Balance is judged by standardized mean differences (treated-group SD), not p-values.

**Report the patients the caliper dropped** (18 of 751 here): they are the sickest treated patients, and
dropping them moves the target towards the ATE.

## Python

No Python package documents MatchIt-compatible defaults, so the greedy algorithm is written out. It
reproduces MatchIt's matched set exactly on this fixture, and the control file replays the greedy order
by brute force.

## Verification

| Engine | Status |
|---|---|
| R (MatchIt) | executed in CI |
| Python | executed in CI |
