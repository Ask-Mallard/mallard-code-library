# Bootstrap percentile and BCa confidence intervals

The ratio of mean length of stay between two hospitals, from skewed data: a statistic with no simple
standard error, which is what the bootstrap is for.

| This fixture (150 patients per hospital) | Value |
|---|---|
| true ratio of means | 1.65 |
| estimate | 1.68 |
| percentile 95% CI | 1.45 to 1.94 |
| **BCa 95% CI** | **1.45 to 1.94** (almost the same here: bias and acceleration are near 0) |

## How

1. **Resample as the data were sampled**: patients within each hospital (a stratified bootstrap),
   2,000 times.
2. **Percentile interval**: the 2.5th and 97.5th percentiles of the bootstrap values.
3. **BCa** (Efron 1987): shifts those percentiles by a bias correction (the share of bootstrap values
   below the estimate) and an acceleration (from the jackknife). Report it by default. Its jackknife
   acceleration needs a smooth statistic: for a median it is exactly 0, which is why this entry uses a
   ratio of means.

Over 100 simulated datasets both intervals covered the true ratio 90 to 91% of the time: skewed means
converge slowly at this sample size, and BCa does not fix that.

## Reproducible across languages

The resamples come from a shared "minimal standard" generator (Park-Miller, multiplier 48271), written
identically in both files, so R and Python agree exactly; the control file checks it against the C++
standard's `minstd_rand` check value. In practice use `boot::boot` or `scipy.stats.bootstrap` with a
seed.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
