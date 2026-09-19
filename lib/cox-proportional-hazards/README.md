# Cox proportional hazards regression

Estimates a hazard ratio while conditioning the baseline hazard away, so no shape is assumed for
it. The proportional hazards assumption is what buys that, and it is an **assumption**: the R file
checks it (`cox.zph`) rather than asserting it.

## The default this entry exists to pin

**R and lifelines both use Efron.** Breslow, the classic alternative, disagrees measurably whenever
event times tie, and is biased toward the null under heavy ties.

Follow-up recorded in whole days or months ties constantly, so this is not an edge case. Every file
in this entry names its tie method explicitly, including R, where Efron is already the default,
because a default that is not written down is one package upgrade from changing the answer.

## The fixture, and a correction worth recording

Times are drawn from an exponential model whose hazard is `baseline * exp(Xb)`, which is
proportional hazards by construction, so the true log hazard ratio is exactly `log(1.6)`.

**The first version of this fixture could not test the thing the entry exists to pin.** It drew
continuous times and produced 714 tied times, every one of them a *censoring* tie at the horizon.
Tied censoring times never invoke a tie-handling method; only tied *event* times do. So Efron and
Breslow would have agreed exactly and the pin would have been decoration.

Times are now rounded to a 0.1 grid, the way real follow-up is recorded. **664 of 785 events now
share a time with another**, and the tie method genuinely matters.

Sanity check needing no model: the crude rate ratio (events per person-time) is 1.5132 against a
true hazard ratio of 1.6000. Lower than the truth, as expected, since a crude rate ratio ignores
the covariate.

## Competing events are out of scope, deliberately

This entry assumes no competing event. Where death competes with the event of interest, a
cause-specific or Fine-Gray model is required, and a Kaplan-Meier curve is not merely off-message:
it censors the competing death and overstates cumulative incidence. `competingEvents` is `false` in
`meta.json` so a plan typing it `true` is never routed here.

## Verification

| Engine | Status |
|---|---|
| R (`survival`) | executed in CI |
| Python (`lifelines`) | executed in CI |
