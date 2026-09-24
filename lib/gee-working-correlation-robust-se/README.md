# Marginal logistic regression by GEE on clustered binary data

Estimates a **population-average** odds ratio when patients within a clinic are more alike than
patients across clinics.

## This is not the same number as the mixed-effects entry

`mixed-effects-logistic-clustered` estimates the **cluster-conditional** odds ratio: the effect for
a patient, holding their clinic fixed. GEE estimates the **marginal** one: the effect in the
population as a whole. On identical data those are different numbers, and the marginal one is
attenuated toward the null.

Neither is a better version of the other, and a plan has to say which it wants. That is also why
this fixture generates the marginal model **directly** rather than through a random intercept: a
random-intercept generator has no closed-form marginal odds ratio, so an entry built that way would
be asking GEE to recover a parameter it does not estimate and reporting the attenuation as error.

## The defaults this entry exists to pin

### The working correlation

| Language | Default working correlation |
|---|---|
| R `geeglm` | **independence** |
| Python `statsmodels.GEE` | **independence** |

Both default to independence, and no output names the structure it used. Every file here writes
`exchangeable` out explicitly, because a default that is not written down is one release from
changing — and the wrong structure is not a rounding difference: it changes both the estimate and
the standard error, as the numbers below show.

### The standard error, and why it interacts with the working correlation

The obvious code in both R (`summary(geeglm)`) and Python (`.bse`) reports the **sandwich**
(empirical) standard error. Its alternative, the **model-based** error, is the **smaller** one, so
reporting it by mistake makes a result look more certain than it is rather than less — the same
shape as the overdispersion trap in `poisson-rate-regression-overdispersed`: **confidently wrong
rather than visibly wrong.**

**How much smaller depends on the working correlation, and that is why this entry reports four
standard errors rather than two.** A model-based variance computed under an EXCHANGEABLE structure
already carries the within-clinic correlation and lands near the sandwich. One computed under
INDEPENDENCE does not carry it at all.

So the two defaults interact, and the severe case is a pair **no single package defaults to**.
Measured in CI on the committed fixture, `geepack`:

| working correlation | model-based SE | sandwich SE | ratio |
|---|---|---|---|
| exchangeable | 0.1188 | 0.1315 | **1.11** |
| independence | 0.0944 | 0.1685 | **1.79** |

The mild row is an exchangeable structure with a model-based variance. R and statsmodels default to
independence with a sandwich, which is also fine. The severe row is independence with a model-based
error, and no package defaults to that pair: it is what a careless edit produces by taking the
model-based error under the default independence structure, and what nobody would write on purpose.

All four are separate harness keys, so an engine reporting one as another fails the agreement check
rather than passing with a plausible number.

**Read the table the right way round.** The mild row is not reassurance: a model-based variance
computed under the *right* structure is close to the sandwich because it has already absorbed the
correlation. It is the pair — independence structure, model-based error — that is 44% too narrow,
and no package produces that pair on its own. A careless edit that pins the model-based error under
the default independence structure does.

### And a silent one in R that has nothing to do with statistics

**`geeglm` requires the rows of a cluster to be contiguous.** It does not check. Unsorted input is
treated as many small clusters, which collapses the sandwich toward the model-based error — the
exact failure the file is written to avoid, arriving through the back door. `r.R` sorts
unconditionally even though the committed fixture already arrives sorted. `statsmodels` does not
need this, which is itself a difference a port from R would lose.

## What the fixture is built to make testable

**Exposure has a clinic-level component**, and this is the choice the entry turns on. Measured on a
working-independence fit — the structure both languages default to — the
sandwich-to-model-based ratio for the exposure coefficient is **1.022** at a clinic-level exposure
spread of 0.9: the demonstration would have been a rounding difference. At 2.5 it is 1.785 on the
committed seed and above 1.3 on 98% of 200 calibration seeds. Some clinics prescribe far more than
others, which is both realistic and what puts the exposure contrast between clusters where the
clustering can reach it.

**Cluster sizes vary from 8 to 40.** The two working correlations coincide exactly in one classic
case — equal cluster sizes with a purely **cluster-level** covariate — and this fixture is away from
it on both counts, so the structure reaches the estimate and not only the standard error. Measured:
`geepack` gives **0.5874** exchangeable against **0.5727** independence.

That is deliberately the weaker phrasing. Unequal cluster sizes do not *prove* the two estimates
must differ; they remove the condition under which the two must agree. The number above is what
establishes that they do differ here.

Size is drawn independently of exposure and outcome. Cluster size that **carries information** about
the outcome makes the two working correlations target genuinely different estimands — real, and not
this entry's subject.

**Ninety clusters**, not twelve. The sandwich is a large-sample estimator and is biased downward
when clusters are few; fewer than about 30 clusters calls for a named small-sample correction, and
fewer than 15 makes one essential. An entry demonstrating the sandwich on a dozen clusters would be modelling the practice
the product warns about. Mancl-DeRouen, Kauermann-Carroll and Fay-Graubard corrections are not in
`geepack` and belong to a separate entry.

The latent within-clinic correlation is 0.25; the resulting **binary** within-clinic ICC is 0.187.
Those are two different numbers and the fixture prints both, because a reader looking for "the ICC"
will find the second and the generator is written in terms of the first.

## The two claims

**Agreement** is between `geepack` and `statsmodels` on the same rows: coefficients, both standard
errors under both working correlations, and the estimated working correlation itself.

`alpha_exchangeable` was included **deliberately** as the key most likely to differ — a check whose
keys are chosen to avoid the awkward ones is calibrated against nothing — **and it does differ.**

| key | `geepack` | `statsmodels` | spread |
|---|---|---|---|
| `alpha_exchangeable` | 0.1849411436 | 0.1847036384 | 2.375e-04 |
| `exposure_se_naive` | 0.1187875776 | 0.1172338035 | 1.554e-03 |

Neither is a defect. The exchangeable α is a **moment estimator that is not uniquely defined**, and
two packages are entitled to different residual scalings; the model-based variance is computed *from*
α, so it inherits the difference.

**What makes that argument checkable rather than a convenient story** is where the disagreement is
*not*. The two engines agree to ten decimals on every key of the independence fit — 0.5726722441,
0.1685252182, 0.0944282807 — and on the exchangeable fit's coefficients and sandwich error. An
independence fit has no α. The disagreement appears exactly where α enters and nowhere else.

So those two keys carry their own agreement tolerances (1e-3 and 5e-3), argued in `expected.json`
rather than quietly widened. They are still 185x and 5x below the gap a **wrong** working
correlation produces (α = 0 against 0.185; a naive SE of 0.0944 against 0.1188), so they admit two
moment estimators and reject a structurally different fit.

**Recovery** is against the fixture's own marginal log odds ratio, `log(2) = 0.6931`, which is exact
by construction rather than a recording of what this code produced. Tolerance 0.50, measured over
200 alternative seeds fitted with a stdlib IRLS logistic regression that is not one of the library's
engines: miss median 0.1190, 95th percentile 0.3283, maximum 0.3978. The committed seed misses by
0.1205.

**That claim is weak and the entry says so.** 0.50 is most of the effect being estimated. With
exposure allocated largely at the clinic level the effective sample size is closer to 90 clinics
than 2169 patients, and the sampling distribution really is that wide — narrowing the tolerance
would fail correct implementations. What carries this entry is agreement, plus a
sandwich-to-model-based ratio that is a property of the design rather than of the seed.

## What is not in this entry

**A small-sample correction.** Not in `geepack`, and a separate entry.

**A conditional odds ratio.** That is `mixed-effects-logistic-clustered`, and it is a different
parameter rather than a different way of getting this one.
