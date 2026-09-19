# Rate regression on overdispersed counts with unequal follow-up

Estimates an **incidence rate ratio**: events per person-time, not events per patient. With
follow-up ranging from six months to five years those are different questions, and only the offset
makes it the first one.

## Two defaults this entry pins

### 1. The offset

| Language | What it takes |
|---|---|
| R | `offset(log(person_time))` — already logged, and **inside the formula**, which fixes the coefficient at 1 |
| Python | `offset=np.log(person_time)` — already logged |

The offset must be the **log** of the person-time. Passing the person-time itself — which some tools'
`exposure()`-style keywords take, logging it for you — fits a model with `exp(person_time)` where
`person_time` belongs. **It runs, it converges, and it answers a different question**, with nothing in
the output to say so.

In R the equivalent slip is writing `+ log(person_time)` as a covariate instead of `offset(...)`.
That estimates the coefficient rather than fixing it at 1, which is a different model — sometimes a
defensible one, never the one you meant when you wrote "rate".

**Follow-up varies here for exactly this reason.** With equal person-time the offset can be dropped
and the rate ratio comes out right anyway, so a fixture with constant follow-up cannot test the
thing most easily got wrong.

### 2. The standard errors, because the counts are overdispersed

The fixture's counts have a **variance 2.4 times their mean**. Under that, Poisson regression's
coefficients are **still consistent** and its model-based standard errors are materially too small.

That is the failure that matters: the estimate is right and the interval is too narrow, so the
analysis is confidently wrong rather than visibly wrong. Both standard errors are reported as
separate harness keys, so an engine reporting the model-based error as robust differs on
`exposure_se` while agreeing on `exposure_log_rr` — which is precisely what the two keys exist to
separate.

An equidispersed fixture could not test any of this: under a true Poisson the naive and robust
errors converge and the pin would be decoration.

## theta and alpha are reciprocals, and every package picks one

| Package | Reports | Variance parameterised as |
|---|---|---|
| R `MASS::glm.nb` | **theta** | mu + mu²/theta |
| statsmodels `nb2` | **alpha** | mu + alpha·mu² |

The fixture's gamma variance is 0.5, so alpha is 0.5 and theta is 2. **Reading R's theta as if it
were alpha describes data four times less dispersed than it is.** Every file here reports alpha; the
R file converts rather than emitting a differently scaled number under the same name.

## The negative binomial keys are the ones to watch

`nb_exposure_log_rr` and `nb_alpha` come from two independent optimisers maximising the same
likelihood, not from a closed form. This library has already measured a pair like that: `clogit` and
`ConditionalLogit` agreed to 5.6e-6, which is why the tolerance here is 1e-4 and not 1e-6. The
dispersion is the parameter most likely to sit at the edge of it.

It is compared anyway. A cross-language check whose keys are chosen to avoid the awkward ones is
calibrated against nothing, and a disagreement here would be worth knowing rather than worth
avoiding.

## The recovery tolerance, and what it is an approximation to

Measured over 300 alternative seeds using the **crude** log rate ratio — total events over total
person-time in each arm, which needs no model and no scientific stack. Median miss 0.0315, 95th
percentile 0.0872, maximum 0.1651; the tolerance is 0.2 and none of the 300 exceeded it. The
committed fixture's crude miss is 0.0391.

The crude estimator approximates the adjusted one's variability rather than reproducing it. The log
link is collapsible and the covariate is independent of exposure, so both estimate the same
parameter, and adjustment can only reduce the variance — so the tolerance is if anything looser than
the adjusted estimator needs. That is stated rather than left to be assumed.

## What this entry does not do

**It is not cluster-robust.** The sandwich here corrects for the variance not matching the mean and
still assumes patients are independent. Where patients share a ward or clinic, or where one patient
contributes several rows, the correction has to be by cluster.

**It does not choose between Poisson-with-robust-errors and the negative binomial.** Both are fitted
and both are defensible: the first makes no assumption about the shape of the extra variation, the
second models it and gives a fully specified likelihood. They answer the same question and disagree
about how much you are willing to assume.

**Zero-inflation is not addressed.** Overdispersion and an excess of zeros are different problems
with different remedies, and a negative binomial fixes only the first.
