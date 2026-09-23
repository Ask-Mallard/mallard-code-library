# E-value and a simple quantitative bias analysis for unmeasured confounding

## The E-value

The minimum strength of association, on the risk-ratio scale, that an unmeasured confounder would need
with **both** the exposure and the outcome to explain away an observed risk ratio (VanderWeele and Ding
2017): RR + √(RR (RR − 1)) for RR > 1. Report it for the estimate and for the confidence limit closer to 1.
It summarizes robustness; it corrects nothing and says nothing about whether such a confounder exists.

## A quantitative bias analysis

When the confounder can be named and its associations taken from external evidence (its prevalence among
the exposed, p1, and unexposed, p0, and its risk ratio with the outcome, RR_UY), the observed RR is
divided by the bias factor (p1 (RR_UY − 1) + 1) / (p0 (RR_UY − 1) + 1). It assumes no
confounder-by-exposure interaction. Report it over a range of bias parameters, or probabilistically.

| This fixture | RR |
|---|---|
| truth (causal) | 1.50 |
| observed (crude) | 2.10 (95% CI 1.85 to 2.37) |
| E-value, estimate / CI limit | 3.61 / 3.11 |
| **bias-adjusted, with the true bias parameters** | **1.68** |

## Verification

| Engine | Status |
|---|---|
| R (EValue 4.1.4) | executed in CI |
| Python (formula written out) | executed in CI |

The control file reproduces the published worked example: RR 3.9 with lower limit 1.8 gives E-values 7.26
and 3.0.
