# Cox regression with a time-varying exposure

When exposure starts during follow-up (a drug started later, a procedure, a diagnosis), each person's
time is split at the start and the exposure is whatever it was during each interval: the
**counting-process** or (start, stop] form. The hazard ratio then compares **currently** exposed with
currently unexposed person-time.

## Immortal time bias

Coding "ever exposed" at baseline credits the time before the exposure started to the exposed group:
time during which they could not have had an exposed event, because they had to survive event-free to
start. It makes the exposure look protective whatever its true effect.

| Coding | Log hazard ratio on this fixture | Truth |
|---|---|---|
| **current exposure, (start, stop]** | −0.508 | −0.511 |
| ever exposed at baseline | **−1.606** | −0.511 |

## The defaults this entry pins

| Language | Pinned |
|---|---|
| R | `coxph(Surv(start, stop, event) ~ exposed, ties = "efron")` |
| Python | `CoxTimeVaryingFitter` with `start_col` and `stop_col` (Efron ties) |

## The fixture

1500 people followed up to 5 years; exposure starts at an exponential time (rate 0.3 a year); event
hazard 0.15 a year while unexposed, × 0.6 after starting.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
