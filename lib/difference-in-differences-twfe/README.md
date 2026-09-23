# Difference-in-differences with two-way fixed effects (single adoption date)

Some hospitals adopt a policy at one time; others never do. The DiD compares the adopters' change with
the non-adopters' change, which removes both the baseline difference between groups and the trend
common to both.

| This fixture | Estimate | Truth |
|---|---|---|
| **DiD (TWFE)** | **−1.13** (SE 0.17) | −1.2 |
| adopters, before vs after | +0.42 | |
| groups, after period | +1.21 | |

## What it assumes

**Parallel trends**: without the policy, the two groups would have changed by the same amount. Plot
the pre-period means of both groups, or fit an event study, before relying on it.

## Scope: one adoption date

With a single adoption date and a balanced panel, the TWFE coefficient equals the 2×2 DiD of group
means exactly (the control file checks it). With **staggered adoption** and effects that change over
time, TWFE compares late adopters with early ones and can be badly biased, even to the wrong sign
(Goodman-Bacon 2021). Use a heterogeneity-robust estimator there (Callaway and Sant'Anna; R `did`).

## The defaults this entry pins

Hospital and quarter fixed effects; SE clustered on hospital (HC1). With fewer than about 30 clusters
see `cluster-robust-cr2-small-sample`.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
