# ANCOVA: a trial effect adjusted for the baseline value

A randomized trial measures the outcome at baseline and at follow-up. Three analyses are unbiased,
because allocation is random, and they differ in precision:

| Analysis | Assumes the baseline carries forward with coefficient | SE on this fixture |
|---|---|---|
| follow-up only | 0 | 0.86 |
| change score (follow-up − baseline) | 1 | 0.71 |
| **ANCOVA** (follow-up on treatment and baseline) | estimated (0.62 here) | **0.61** |

ANCOVA is the most precise because it estimates the carry-over instead of assuming it.

## Rules

- **Adjust only for what was measured before randomization.** A variable measured after treatment
  starts can be affected by it.
- Model-based standard errors suit equal allocation with similar spread in both arms. With clearly
  unequal variances or unequal allocation, use robust errors (`linear-regression-robust-se`).

## The fixture

200 per arm; baseline ~ N(50, 10); follow-up = 20 + 0.6 × baseline − 3 × treated + noise (SD 6).
True effect −3. Observed −3.90.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
