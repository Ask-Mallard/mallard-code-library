# Bivariate random-effects meta-analysis of diagnostic accuracy (Reitsma)

Sensitivity and specificity vary together across studies (a lower positivity threshold raises one and
lowers the other), so they are pooled **jointly**: the bivariate model treats each study's logit
sensitivity and logit specificity as bivariate normal (Reitsma et al. 2005), fitted by REML.

| This fixture (20 studies) | Truth | Estimate |
|---|---|---|
| summary sensitivity | 0.818 | 0.808 |
| summary specificity | 0.881 | 0.879 |
| between-study SD (logit sensitivity / specificity) | 0.50 / 0.60 | 0.62 / 0.53 |
| correlation | −0.4 | −0.02 |

Report the summary point with its confidence region and a prediction region (or summary ROC curve):
between-study variation in accuracy is usually large and belongs in the conclusion. With 20 studies the
correlation is poorly estimated.

## R only

Python's standard stack has no bivariate random-effects meta-analysis. The control file refits the same
model with `metafor::rma.mv` and requires `mada::reitsma` to match it.

## Verification

| Engine | Status |
|---|---|
| R (mada 0.5.12) | executed in CI |
| Python | not applicable (see above) |
