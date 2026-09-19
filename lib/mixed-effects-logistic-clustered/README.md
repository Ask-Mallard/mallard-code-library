# Mixed-effects logistic regression for clustered binary outcomes

Estimates a **cluster-conditional** odds ratio with a random intercept per cluster.

## There is no `python.py` here, and that is the answer rather than a gap

The standard Python stack has no frequentist GLMM equivalent to `lme4::glmer`. Three candidates
exist and each was rejected for a specific reason:

| Candidate | Why not |
|---|---|
| `statsmodels` `BinomialBayesMixedGLM` | Variational Bayes: a different estimator targeting a different thing. Agreement with `glmer` would be meaningless either way. |
| `statsmodels` `GEE` | A **population-average** model. It targets an attenuated quantity, not the same number by another route. |
| `pymer4` | Wraps `lme4` through `rpy2`, so it *is* R. Two calls to one engine agreeing is not independent evidence. |

So this entry makes the **recovery** claim only, and `harness/check.py` refuses to state the
agreement claim with fewer than two engines. The weaker guarantee is reported as weaker.

Substituting GEE to make the entry look complete would have been the worst option available: it
would silently swap the estimand, and the agreement check would then compare two different
parameters and report the gap as a defect.

## Conditional is not marginal

A mixed model answers "what is the effect for a given cluster". A GEE answers "what is the effect
averaged over the population". Under a nonlinear link those are **different quantities, not
different spellings**, and the marginal one is attenuated toward 1.

The fixture shows it without needing either model: the crude **marginal** odds ratio is **1.379**
against a true **conditional** odds ratio of **2.000**.

Which one answers your question is a choice to make and state. It is not an inefficiency to be
tolerated, and it is why "mixed-effects versus GEE" is one of the five comparisons in the Learn
plan.

## The fixture

40 clusters of 30, exposure allocated **at the cluster level**, latent-scale ICC 0.163.

Cluster-level allocation is deliberate. With individual-level exposure inside clusters the random
effect is nearly orthogonal to it, and an analysis that ignored clustering entirely would still
land near the truth: the fixture would not test what the entry is for.

## Sizing

`calc: "simulation"`. A clustered design is sized by simulation, not a closed form, and reporting
an unclustered two-group formula or a per-variable floor as *the* method for a clustered design is
a documented error. The alternative to a wrong number is no number.

## The recovery tolerance, and the requirement it turned out could not be met

This entry shipped with a **provisional** tolerance of 0.45, chosen so it would stay narrower than
the gap between the conditional log odds ratio (0.6931) and the crude marginal one (0.3213) — the
reasoning being that a wider tolerance would stop distinguishing the two estimands, which is the one
thing this entry must not do.

**Calibration showed those two requirements cannot both be met at 40 clusters.** Fitted over 300
alternative seeds with a stdlib random-intercept MLE by 30-node Gauss-Hermite quadrature — not
`lme4`, not `statsmodels`, and it reproduces `glmer`'s answer on the committed fixture to five
decimals — the miss distribution is:

| | median | 95th | max |
|---|---|---|---|
| conditional log OR | 0.1927 | 0.5598 | 0.8859 |
| between-cluster SD | 0.0832 | 0.2414 | 0.4033 |

The sampling distribution needs **0.89**; the estimand gap allows **0.37**. At 0.45, roughly one run
in ten of a **correct** implementation would have failed. Narrowing a tolerance below the sampling
distribution does not sharpen a check — it makes it fire on correct work, which is how a guard
earns the reputation that gets it deleted.

**So the estimand distinction is carried by `cluster_sd` instead**, at its own tolerance of 0.45.
A population-average model does not produce a between-cluster SD at all: there is no such parameter
in it. Recovering 0.8 within 0.45 is therefore a claim only a conditional model can make, and it
separates the two estimands **structurally** rather than by a numerical margin the data cannot
support. That is also why the harness grew per-key recovery tolerances: 0.90 applied to a parameter
whose entire miss distribution tops out at 0.40 would pass whatever it was handed.

## Verification

| Engine | Status |
|---|---|
| R (`lme4`) | executed in CI |
| Python | **not applicable** — see above |
