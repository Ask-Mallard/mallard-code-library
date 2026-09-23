# Comparing two diagnostic tests on the same patients

When both tests are done on everyone, their results are correlated and the comparison must be paired.

| This fixture | Truth | Estimate |
|---|---|---|
| AUC difference (A − B) | 0.112 | 0.102 (paired SE 0.012; unpaired 0.017) |
| sensitivity difference at threshold 0.6 | 0.186 | 0.152 (McNemar χ² 42.6) |

## The defaults this entry pins

- **AUCs**: DeLong's test for two correlated ROC curves (`roc.test(paired = TRUE)`), using the
  covariance of the two AUCs. Ignoring it overstates the SE when the markers are positively correlated.
- **At a threshold**: sensitivities compared among patients with disease (specificities among those
  without) by **McNemar's test** on the discordant pairs, without continuity correction. A two-proportion
  test ignores the pairing. The threshold must be **prespecified**.
- Direction pinned for both markers: higher values mean disease.

## Verification

| Engine | Status |
|---|---|
| R (pROC 1.19.1) | executed in CI |
| Python (paired DeLong written out) | executed in CI |

Python's paired z equals pROC's `roc.test` statistic.
