# Small-sample cluster-robust inference: CR2 with Satterthwaite df

A patient-level regression in a trial with **few clusters of unequal size**. The ordinary cluster-robust
(sandwich) standard error is too small with few clusters, and a normal reference makes the interval
narrower still. **CR2**, the Bell-McCaffrey bias-reduced sandwich, with **Satterthwaite degrees of
freedom**, keeps coverage close to nominal (Pustejovsky and Tipton 2018).

| Standard error of the treatment effect | This fixture |
|---|---|
| naive OLS (patients independent) | 0.55 |
| CR1 (Stata-style cluster sandwich) | 1.27 |
| **CR2** | **1.33**, with **8.98** df (not G − 1 = 11) |

## The defaults this entry pins

R: `clubSandwich::coef_test(vcov = "CR2", test = "Satterthwaite")`. Python has no package for it
(statsmodels' cluster option is CR1-style with G − 1 df), so `python.py` writes CR2 and the
Bell-McCaffrey df out. It matches clubSandwich to 1e-10, including the fractional df.

## The fixture

12 clinics (6 per arm) of 20 to 60 patients; outcome 20 + 4 × treated + 0.1 × (age − 50) + clinic
effect (SD 1.5) + noise (SD 5). True effect 4; observed 3.18.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
