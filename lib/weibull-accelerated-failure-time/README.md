# Weibull regression: time ratio and hazard ratio

A parametric survival model assumes a shape for the hazard: for a Weibull, steadily rising (shape > 1)
or falling (shape < 1). In return it gives a **time ratio** ("treatment stretches survival times by
38%") and can extrapolate beyond follow-up, which Cox regression cannot. The extrapolation is only as
good as the assumed shape.

## Time ratio, not hazard ratio

The model's coefficient is on the log **time** scale: positive means longer survival. For a Weibull,
the equivalent log hazard ratio is −coefficient × shape, with the opposite sign.

| Quantity | This fixture | Truth |
|---|---|---|
| log time ratio (the coefficient) | +0.325 | +0.4 |
| shape | 1.495 | 1.5 |
| log hazard ratio (−coefficient × shape) | −0.487 | −0.6 |

Reading the coefficient as a log hazard ratio gets both the sign and the size wrong.

## Parametrizations

R's `survreg` reports `scale` = 1 / shape. lifelines' `WeibullAFTFitter` models the shape (`rho_`) on
the log scale. Both files convert to the shape explicitly.

## Agreement

lifelines stops about 5e-6 short of the maximum likelihood here, and its `fit_options` did not tighten
that. R's survreg is at the maximum, confirmed by an independent optimisation, so the entry uses the
library's standard 1e-4 agreement tolerance. The difference is far below anything clinical.

## The fixture

600 people; Weibull shape 1.5, scale 5 years (control), treatment multiplies times by e^0.4; follow-up
ends between 2 and 8 years.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
