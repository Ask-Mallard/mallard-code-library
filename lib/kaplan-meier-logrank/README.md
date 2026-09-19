# Kaplan-Meier estimation and the log-rank test

Estimates survival curves from partially observed follow-up, and compares two of them.

**Neither of these produces an effect estimate.** Kaplan-Meier gives a curve and the log-rank gives
a test; a hazard ratio is the Cox entry. A plan that needs a number to report needs one of those,
and a log-rank p quoted beside no estimate is the shape this entry must not encourage.

## The default this entry exists to pin

**R and Python put a different default band around the same curve, and neither says so in its
output.**

| Language | Default confidence transformation |
|---|---|
| R | `survfit` → **`"log"`** |
| Python | lifelines → exponential Greenwood, which **is** the log-log transform |

Every file here pins **log-log**, because it is the one that cannot leave `[0, 1]`. The plain log
transform produces upper limits above 1 in the tail, where a survival probability cannot go. R is
the language that has to be told, because its default is the plain log; a default that is not written
down is one release from changing.

**And the test has its own version of the problem.** `survdiff`'s `rho = 0` is the log-rank; `rho = 1`
is Peto-Peto, which weights early differences more heavily, so every file names `rho = 0` explicitly.

## What the log-rank is for, and what it is not

It has most power against **proportional hazards**, which is what this fixture generates — it is the
score test of exactly the model the Cox entry fits.

It is **not** a general "are these two curves different" test. Where curves **cross** it can return a
p near 1 while the arms are plainly different, because early and late differences cancel. That is a
property of the test rather than a fault, and it is why every file here plots or lists the curves
beside the p value.

## Three things the fixture is built to make testable

**Censoring is heavy and independent** — 27% of the 900 patients, from loss to follow-up as well as
the horizon. Kaplan-Meier exists to use partial follow-up; a fixture where nearly everyone is
observed to fail would be estimated about as well by ignoring censoring, and the estimator would not
be doing any work.

**Event times tie** — 607 of 656 events share a time with another, because follow-up is recorded on
a half-month grid, as real follow-up is. Ties are what make the step heights depend on a consistent
convention.

**The median exists in both arms**, and the generator checks that before writing the file. A
Kaplan-Meier median is undefined when the curve never reaches 0.5, and every language spells that
differently — R gives `NA`, lifelines gives `inf`. An entry whose fixture hit that case would be
comparing two spellings of "undefined" and calling it agreement.

## The two claims

**Agreement should be exact, not approximate**, for the medians and the step values: both engines
evaluate the same step function on the same rows, so any difference is a convention difference
rather than arithmetic.

**The interval keys are the ones to watch.** This is the first entry in the library to compare an
*interval* rather than a point estimate — R's `conf.type="log-log"` against lifelines' exponential
Greenwood. Those are the same transform, but they are two independent implementations of it. Included
deliberately: a check whose keys are chosen to avoid the awkward ones is calibrated against nothing.

**Recovery is against a closed form.** The generator is exponential, so survival at twelve months is
`exp(-0.06 × 12)` and `exp(-0.06 × 1.8 × 12)` exactly — arithmetic, not a recording of what this code
produced. Tolerance 0.12, measured over 200 alternative seeds: worse-arm miss median 0.0248, 95th
percentile 0.0552, maximum 0.0853. The committed seed's worse miss is 0.0275.

The half-month grid biases the estimate slightly **upward** — an event can move from just before
twelve months to just after — which is a property of recording follow-up on a grid and part of what
the fixture is for.

## What is not in this entry

**Competing events.** Where death competes with the event of interest, a Kaplan-Meier curve is not
merely off-message: it censors the competing death and **overstates cumulative incidence**. The
cumulative incidence function is the estimator, and `competingEvents` is `false` in `meta.json` so a
plan typing it true is not routed here.

**Adjustment.** These are crude curves. Adjusted survival needs a model.

**A hazard ratio.** That is `cox-proportional-hazards`, whose fixture is built the same way for the
same reasons.
