# Sharp regression discontinuity: local linear with robust bias-corrected inference

Treatment is assigned by a threshold on a score (at or above the cut-off: treated). The effect is the
**jump** in the outcome at the cut-off, for patients near it.

| This fixture | Value |
|---|---|
| truth | 1.5 |
| **local linear estimate** | **1.59**, robust 95% CI 1.20 to 2.06 |
| bandwidth h (MSE-optimal) | 0.28 |
| global straight line on each side | 0.66 |

## The defaults this entry pins

`rdrobust` in both languages (same authors): local linear (p = 1) with a triangular kernel within an
MSE-optimal bandwidth, a local quadratic (q = 2) bias correction, nearest-neighbour variance with 3
matches. **Report the conventional estimate with the robust bias-corrected interval**: at the
MSE-optimal bandwidth the conventional interval is too narrow because the estimate still carries
smoothing bias (Calonico, Cattaneo and Titiunik 2014).

Global polynomials are not a substitute: a straight line on each side misses badly when the outcome
curves, and high-order polynomials are unstable near the edge (Gelman and Imbens 2019).

## Checks to report

Density of the score around the cut-off (manipulation: `rddensity`), and covariate balance at the
cut-off. This fixture has neither problem.

## Verification

| Engine | Status |
|---|---|
| R (rdrobust 4.0.0) | executed in CI |
| Python (rdrobust 2.0.0) | executed in CI |

The control file refits the local linear regression by hand at rdrobust's bandwidth.
