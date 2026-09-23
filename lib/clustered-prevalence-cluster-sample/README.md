# Prevalence from a cluster sample, with a design-based interval

Clinics, wards or practices were sampled, and then every eligible patient in each. Patients in one
clinic resemble each other, so the sample holds less information than its patient count suggests.
A Wilson interval over the patients treats them as independent and is too narrow.

## Two different "prevalences"

| Quantity | Answers | This fixture |
|---|---|---|
| **Patient-level prevalence** (total cases / total patients) | what share of patients have the condition | 0.15 |
| Average of the clinic prevalences | what share a typical clinic has | 0.20 |

They are the same only when clinic size is unrelated to prevalence. This entry estimates the first,
which is almost always the question. The fixture makes small clinics run high and large clinics run
low, so averaging the clinic percentages is caught by the recovery check on the committed rows.

## The defaults this entry pins

| Language | Default | Consequence of relying on it |
|---|---|---|
| R `svydesign` | `ids = ~1` if the clusters are not named | patients treated as independent; SE too small |
| R `confint` on a survey estimate | `df = Inf` | a normal interval, too narrow with few clusters |
| Python | no survey package in the standard stack | the linearization is written out, with the `m/(m-1)` factor and a t reference on clusters − 1 |

## Scope

One-stage cluster sample, clusters drawn with equal probability, with-replacement variance, no
finite population correction, complete data. Stratified designs use
`complex-survey-domain-prevalence`. With fewer than about 30 clusters, the t reference helps, but the
plan should still say the interval rests on few clusters.

## The fixture

40 clinics from 400, sizes 10, 30 or 60 with prevalence about 30%, 20% and 10%. 1,270 patients.
Observed patient-level prevalence 0.1795 against a true 0.15, about 2.1 standard errors. On these
rows the SE that ignores clustering is 1.28 times smaller than the design-based one.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |

`harness/descriptive_examples.test.py` checks the SE against the variance computed from cluster
counts and measured the recovery tolerance over 100 seeds.
