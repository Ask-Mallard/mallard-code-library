# Binary prediction: discrimination, calibration and Brier score

Evaluates fixed-model validation predictions in independent participants. Unweighted complete-data point estimates, with no fitting or tuning of the prediction model on the validation sample.

The prediction model is fixed before validation outcomes are examined. No updating or threshold selection is performed here.

Calibration-in-the-large uses log odds as an offset, fixing its slope to one. The joint calibration model estimates a separate intercept and slope.

AUC direction is fixed: higher predictions mean outcome=1. Tied predictions receive half credit.

Predictions must be strictly between zero and one; the example rejects invalid values rather than silently clipping them.

These are point estimates for uncensored, independent, unweighted observations; real evaluation also needs uncertainty, a calibration curve and a clinical-use assessment.

R and Python are executed on the same synthetic fixture. Agreement and recovery are separate checks. Generated adaptations and analyses on user data have not been executed by these checks.

Sources: https://www.bmj.com/content/384/bmj-2023-074819 and https://www.bmj.com/content/384/bmj-2023-074820
