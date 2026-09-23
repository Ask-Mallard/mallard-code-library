# Calibration of a prediction model: slope, calibration-in-the-large and O/E

Discrimination (AUC) is not calibration: a model can rank patients well and still predict 30% for
patients whose risk is 15%. Validating a published model in new patients means checking both.

| This fixture | Truth | Estimate |
|---|---|---|
| calibration slope | 0.667 | 0.676 |
| calibration-in-the-large | −0.727 | −0.767 |
| observed / expected | 0.670 | 0.654 |
| AUC (for contrast) | | 0.734 |

## Three numbers, each from its own model (Van Calster et al. 2019)

- **Calibration slope**: the coefficient of the linear predictor in `glm(y ~ lp)`. Below 1 means
  predictions too extreme, the signature of overfitting.
- **Calibration-in-the-large**: the intercept of `glm(y ~ offset(lp))`, with the slope **fixed at 1**.
  The intercept of the slope model is a different quantity (−0.86 here).
- **Observed/expected**: total events over the sum of predicted risks.

Plot the **calibration curve** as well (loess or a spline). It is not compared across languages here,
because R's loess and Python's lowess are different smoothers. The Hosmer-Lemeshow test is not a
substitute.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
