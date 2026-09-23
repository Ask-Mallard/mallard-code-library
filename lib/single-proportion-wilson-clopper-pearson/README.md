# A single proportion with Wilson and Clopper-Pearson intervals

A prevalence from a simple random sample: the share of people with a condition, with a 95%
interval. This is the one case the Wilson and Clopper-Pearson intervals are for.

**That licence does not extend anywhere else.** An odds ratio, a risk ratio or an adjusted model
estimate takes a profile-likelihood or Wald interval. Patients sampled in clinics, wards or practices
are not independent, and their prevalence takes a design-based interval
(`clustered-prevalence-cluster-sample`).

## Which interval

| Interval | When |
|---|---|
| Wilson | the usual report; close to nominal coverage even near 0 and 1 |
| Clopper-Pearson ("exact") | when a conservative interval is wanted, for example with very few events; always at least 95% coverage, so wider |
| Wald (p ± 1.96 SE) | not used: too narrow near 0 and 1, and can run below 0 |

## The defaults this entry pins

| Language | Default | Consequence of relying on it |
|---|---|---|
| R `binom.confint` | `methods = "all"` | eleven intervals; you report whichever row you indexed |
| Python `proportion_confint` | `"normal"` | the Wald interval |
| Python `proportion_confint` | no method named `"exact"` | Clopper-Pearson is `method="beta"` |

## The fixture

400 people, true proportion 0.12. The proportion is kept away from 0.5 on purpose: near 0.5 the
three intervals nearly coincide and a silent fall-back to Wald would pass. Here the Wald lower bound
differs from Wilson's by 0.0035, which the controls assert. Observed: 53 events, 0.1325.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |

`harness/descriptive_examples.test.py` also checks the Python output against the closed-form
Wilson interval and the beta-quantile Clopper-Pearson interval, and measured the recovery tolerance
over 100 seeds.
