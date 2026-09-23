# Self-controlled case series (conditional Poisson)

Only cases are sampled. The event rate in a **risk window** after exposure is compared with the rate at
other times in the **same person's** observation period, so everything fixed over that period (sex,
genes, underlying frailty) is removed, measured or not. It suits acute events after a transient exposure,
such as a vaccine.

| This fixture | log IRR | IRR |
|---|---|---|
| truth | 1.099 | 3.0 |
| **SCCS with age** | **1.204** (SE 0.096) | 3.33 |
| SCCS without age | 1.662 | 5.27 |

## Age must be in the model

Exposure is usually given at particular ages, and event rates change with age. Here vaccination falls in
the high-risk age band; leaving age out credits the age effect to the vaccine.

## The model

A Poisson regression with **one fixed effect per case** and log(interval length) as the offset. It has
the same estimate and standard error as the SCCS conditional (multinomial) likelihood; the control file
fits that likelihood directly. The R `SCCS` package's `standardsccs()` fits the same model.

## Assumptions

Events do not change the probability of later exposure or the end of observation (death after the event
breaks this); recurrent events are independent, or only the first event is used.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI (`wls_method="qr"`: the default SVD solver failed on 1 of 100 calibration fixtures) |
