# An average treatment effect by inverse-probability-of-treatment weighting

Weights each patient by the inverse of their probability of receiving the treatment they actually
received, so the weighted population looks like one in which treatment was assigned at random.

## The point estimate is the easy part

Every implementation below gets essentially the same number. **Every default standard error is
wrong**, and each one is wrong differently. That is what this entry is for.

### The weight argument means a different thing in each language

| Language | What you write | What the weights are taken to be |
|---|---|---|
| R | `glm(..., family = gaussian, weights = w)` | **prior precision** — correct shape here |
| R | `glm(..., family = binomial, weights = w)` | a **count of trials** — "non-integer #successes", and a standard error as if the study had `sum(w)` participants |
| Python | `WLS(..., weights = w)` / `var_weights` | **variance weights** — correct shape here |
| Python | `GLM(..., freq_weights = w)` | **counts**, same error as R's binomial case |

**R has one argument whose meaning is chosen by the family; statsmodels has two and makes you pick.**
The same intent, spelled more than one way, with a different variance behind each.

### And the sandwich is required, not optional

The model-based error from a weighted fit treats the weights as **known precision**. They are not:
they were estimated from a model fitted to the same data. `sandwich::vcovHC(type = "HC0")` in R and
`cov_type="HC0"` in statsmodels are the defensible choice, and both files report the model-based
error beside the robust one so the gap is visible rather than described.

**It is still conservative, and every file says so.** The sandwich here treats the propensity score
as fixed. Accounting for the fact that it was estimated — an M-estimation stack over both stages, or
a bootstrap of the whole two-stage procedure — gives a **narrower** interval. Conservative is the
safe direction and is not the same thing as right. An entry that reported this interval as exact
would be teaching the error it exists to prevent.

## What the fixture is built to make testable

**The truth is a population quantity, not an estimate.** Both potential outcomes are drawn for every
patient and only the observed one is written out, so

```
ATE = E_L[expit(b0 + b1 + b2·L)] − E_L[expit(b0 + b2·L)] = 0.1457586545
```

is an integral against a standard normal, evaluated by 60-node Gauss-Hermite quadrature to about
1e-12. It exists before any estimator does.

**That choice is load-bearing.** IPTW targets the **ATE**; matching without replacement targets the
effect among the **treated**. On confounded data those are different numbers rather than two routes
to one, so an entry whose truth came from one estimator would be scoring the other for disagreeing.

**Confounding is large.** The crude risk difference is **0.2447** against a true **0.1458** — the
unadjusted analysis overstates the effect by 68%. A fixture where adjustment barely moved the
estimate could not tell correct weighting from none.

**Positivity holds with room, and that is asserted rather than hoped.** The propensity score spans
[0.081, 0.935]; the generator fails if it leaves (0.05, 0.96) or if the largest weight exceeds 15,
so a future edit that walks into the extreme-weight regime breaks rather than passes. Extreme
weights are a real problem in real data, but an entry whose fixture produced them would be testing
each implementation's **trimming** defaults rather than its weighting — a separate decision, and a
separate entry.

The cost of that choice is recorded rather than hidden: a stronger treatment model (`a1 = 0.9`)
gives more confounding, a crude estimate of 0.288, and a largest weight of 23.5. The weaker model
was taken because the entry is about weighting, not about surviving bad weights.

**Three weight numbers, kept apart.** The generator knows the true propensity and prints the weights
it implies (10.5 unstabilized, 5.30 stabilized). The analysis *estimates* the propensity, and its
largest stabilized weight is **4.47** — the harness key `max_stabilized_weight`. Quoting one as the
other would misdescribe the fixture in whichever direction flattered it.

## The two claims

**Agreement** is between R (`glm` + `sandwich`) and `statsmodels` on the same rows: the effect, both
standard errors, the propensity coefficient, and the weight diagnostics.

`ate_se_robust` is the key to watch — R's `vcovHC(type="HC0")` on a weighted `glm` against
statsmodels' `WLS(cov_type="HC0")`. The same estimator mathematically, but two implementations of a
weighted sandwich, and the packages' conventions inside it need not coincide. It is included at the
**ordinary** tolerance deliberately: if they differ, the entry says which and why rather than
loosening the key.

`mean_stabilized_weight` is a **diagnostic, not decoration**. Stabilized weights average to about 1,
and across the 300 calibration seeds this one stays inside 0.998 to 1.002. A run that reports 1.9
has unstabilized weights; one that reports something far from 1 has a propensity model that does not
fit.

**Recovery** is against the fixture's own ATE, tolerance 0.06, measured over 300 alternative seeds
with a stdlib IRLS propensity model and a Hájek weighted risk difference — neither of them one of
the library's engines. Miss: median 0.0094, 95th percentile 0.0312, maximum 0.0443. The committed
seed misses by 0.0107.

**And that tolerance rejects the unadjusted estimator outright.** Across the same 300 seeds the
crude risk difference misses the ATE by **at least 0.0736**, median 0.1185 — every one of the 300 is
outside 0.06. So an implementation that silently dropped the weights, or inverted them, could not
pass this entry's recovery claim. The check separates the analysis from *no* analysis, not merely
from a wild one.

## What is not in this entry

**Trimming, truncation and extreme weights.** The fixture deliberately stays out of that regime.

**Matching.** A different estimand (the ATT) and a different set of traps — chiefly that a naive
standard error on a matched sample ignores the matching.

**Doubly robust estimation.** Augmented IPTW gets one more chance at being right and has its own
failure modes.

**Any claim that the propensity model is correct.** It is correctly specified here *by construction*,
which is exactly what you never have. The entry pins how the weights are used, not how to choose
them — and a propensity model chosen to maximise discrimination can be worse for the estimate than
one that is merely correct, which is why balance and not the c-statistic is the check named in every
file.


