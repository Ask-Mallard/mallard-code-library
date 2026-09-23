# Predictive values and likelihood ratios

| This fixture (prevalence 20%) | Truth | Estimate |
|---|---|---|
| PPV | 0.68 | 0.667 (Wilson 0.618 to 0.712) |
| NPV | 0.96 | 0.960 (0.947 to 0.970) |
| LR+ | 8.5 | 8.03 |
| LR− | 0.167 | 0.168 |
| PPV if prevalence were 5% | 0.309 | 0.297 |

## Predictive values belong to a population

PPV and NPV depend on prevalence: the same test has a PPV of 0.68 here and about 0.30 where the disease
is four times rarer. Carry sensitivity and specificity (or likelihood ratios) to a new setting and
recompute. A case-control sample, whose prevalence is set by design, gives no valid predictive values.

## Intervals match the estimator

PPV and NPV are single proportions: Wilson intervals. Likelihood ratios are ratios of independent
proportions: log-scale intervals (Simel, Samsa and Matchar 1991).

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
