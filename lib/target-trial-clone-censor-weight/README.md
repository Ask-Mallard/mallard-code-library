# Target trial emulation with a grace period: clone, censor, weight

Observational data can answer a trial's question if the analysis emulates the trial: the same eligibility,
the same **time zero** for every strategy, and a comparison of strategies rather than of the treatment
people happened to get.

## The target trial here

- Eligible at interval 0 (time zero), followed for 10 intervals.
- Strategies: **start treatment within 3 intervals** (the grace period) against **never start**.
- Outcome: an event within 10 intervals; contrast: per-protocol risk difference.

## Clone, censor, weight

1. **Clone** every patient into both strategies: at time zero their data are consistent with both.
2. **Censor** each clone when its data stop agreeing with its strategy: the never clone when the patient
   starts, the grace clone at the end of the grace period if the patient has not started.
3. **Weight** the uncensored clones by the inverse probability of remaining uncensored, because starting
   depends on severity, which also drives the outcome.
4. Estimate each arm's risk by the weighted Kaplan-Meier.

| This fixture | Risk difference |
|---|---|
| truth (simulated under each strategy) | −0.214 |
| **clone, censor, weight** | **−0.215** |
| clone and censor, no weights (confounded) | −0.175 |
| ever started against never (immortal time) | −0.249 |

Confidence intervals need a bootstrap of the whole procedure; this entry is checked against its truth only.

## R only

Target trial emulation is included with a truth check only, so this entry makes no agreement claim. The control
file re-implements the procedure in numpy and must reproduce r.R's output, and recomputes the truth with a
second Monte Carlo seed.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | not applicable (see above) |
