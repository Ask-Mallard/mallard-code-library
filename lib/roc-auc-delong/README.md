# ROC analysis: AUC with a DeLong interval and a Youden threshold

The area under the ROC curve (AUC) is the probability that a patient with disease has a higher marker
value than one without. It summarizes discrimination across all thresholds; it says nothing about
calibration or about performance at the threshold a clinic would use.

| This fixture | Value |
|---|---|
| true AUC | 0.802 |
| **estimate** | **0.764** (DeLong 95% CI 0.724 to 0.805) |
| Youden threshold | 0.52 (sensitivity 0.71, specificity 0.70) |

## The defaults this entry pins

- **Direction**: higher values mean disease. pROC's default (`direction = "auto"`) chooses whichever
  direction gives AUC > 0.5, which turns a marker pointing the wrong way into a good one.
- **DeLong** interval on the AUC scale.
- **Youden threshold** among midpoints of consecutive observed values. A threshold chosen on the same
  data is optimistic: over 100 simulated datasets the Youden index there exceeds the true maximum by
  0.015 on average. Prespecify or validate a threshold for use.

## Verification

| Engine | Status |
|---|---|
| R (pROC 1.19.1) | executed in CI |
| Python (DeLong written out) | executed in CI |
