# Random-effects meta-analysis: REML with Hartung-Knapp

Fifteen trials of one treatment, each reporting a 2×2 table. The true effects vary between trials, so
the summary is the **mean of a distribution of true effects**, and how much they vary matters as much as
the mean.

| This fixture | Value |
|---|---|
| true mean log OR (OR) | −0.40 (0.67) |
| **summary OR, Hartung-Knapp 95% CI** | **0.70 (0.59 to 0.82)** |
| τ² (REML), I² | 0.020, 24% |
| 95% prediction interval (OR) | 0.49 to 0.98 |

## The defaults this entry pins

- **REML** for τ² (DerSimonian-Laird underestimates it with few trials).
- **Hartung-Knapp** for the interval of the mean: t on k − 1 df with a rescaled variance.
- A **prediction interval** for the effect in a new trial, always wider than the confidence interval.
- I² is a proportion of variability, not an amount; report τ² too.

Over 100 simulated meta-analyses Hartung-Knapp covered the true mean 92% of the time,
DerSimonian-Laird with a z interval 91%.

## Verification

| Engine | Status |
|---|---|
| R (metafor 5.0-1 in CI) | executed in CI |
| Python (REML written out) | executed in CI |
