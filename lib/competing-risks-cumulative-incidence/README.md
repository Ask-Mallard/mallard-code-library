# Competing risks: the cumulative incidence function

Estimates how many patients **actually have** the event, when a competing event can stop them from
ever having it.

## The error this entry exists to put numbers on

A Kaplan-Meier curve for cause 1 treats a death from the competing cause as **censoring** — as
though that patient could still go on to have cause 1. They cannot. So the naive curve overstates
cumulative incidence, and the overstatement grows with the competing hazard.

On this fixture, where the competing cause is the more common one:

| | 12 months | 24 months | 36 months |
|---|---|---|---|
| Cumulative incidence (truth, closed form) | 0.2257 | 0.3070 | 0.3364 |
| 1 − Kaplan-Meier (the naive curve) | 0.3023 | 0.5132 | 0.6604 |
| overstated by | 0.077 | 0.206 | 0.324 |

**At 36 months the naive curve nearly doubles it.** That is the failure CLAUDE.md's statistical
rules name directly, and this entry is it with numbers on.

**Every file computes the naive curve on purpose**, beside the right one, so the size of the error
is in the output rather than in a comment. It is the line to delete from your own analysis, not to
copy.

## Three spellings, three answers, no warning

This is where the entry earns its keep: in every language the wrong analysis is a *shorter* piece of
code than the right one.

| | Right | Silently counts cause 2 as cause 1 | Silently censors cause 2 |
|---|---|---|---|
| R | `Surv(time, factor(status))` | `Surv(time, status)` | `Surv(time, status == 1)` |
| Python | `AalenJohansenFitter(...event_of_interest=1)` on the full column | — | a 0/1 indicator |

None of those spellings warns you. The R column is the sharpest: three characters of difference
between an estimator that is right, one that pools two different diseases, and one that is
optimistic by a factor of two.

## The fixture has no tied event times, unlike the Kaplan-Meier entry

Deliberate, and for a reason that is a property of one engine rather than a statistical one:
lifelines' `AalenJohansenFitter` does not support tied event times and adds **jitter** when it finds
them, which would make the two engines disagree for a reason that has nothing to do with the
estimator.

This entry's subject is the gap between two **estimators**, so the fixture avoids ties instead of
testing them, and the generator asserts zero ties before writing the file. The Kaplan-Meier entry
does the opposite, deliberately, because there ties *are* the subject.

## The two claims

**Recovery is against a closed form.** With constant cause-specific hazards the cumulative incidence
of cause 1 is `h1/(h1+h2) × (1 − exp(−(h1+h2)t))` — arithmetic, not a recording of what this code
produced. Tolerance 0.06, measured over 120 alternative seeds: worst-of-three miss median 0.0135,
95th percentile 0.0277, maximum 0.0362. The committed seed's worst is 0.0102.

**And that tolerance rejects the naive curve outright**, which is the property worth having: 1 − KM
misses by 0.078, 0.197 and 0.308, so an implementation that quietly censored the competing cause
could not pass this entry's recovery claim.

**The naive value is not in `truth`.** It has a limit — `1 − exp(−h1 t)` — it is simply the wrong
quantity, and putting it under `truth` would read as though the library sanctioned it. It is
compared between engines instead, so a change that quietly fixed or broke the demonstration is still
caught.

## What is not in this entry

**A regression.** Fine-Gray models the **subdistribution** hazard; a cause-specific Cox model
censors the competing event and answers a different question —
aetiology versus prediction. Both are out of scope here and neither is a substitute for the other.

**A comparison between groups.** Gray's test is the analogue of the log-rank on this scale. One
fixture, one curve.

**Sizing.** `calc` is `none`: the catalog has no closed form for a cumulative incidence comparison,
and borrowing the log-rank form would size a different quantity.
