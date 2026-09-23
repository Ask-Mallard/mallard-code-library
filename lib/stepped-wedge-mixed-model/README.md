# Stepped-wedge trial: the Hussey-Hughes mixed model

In a stepped wedge every cluster starts in control and switches to the intervention at a randomized time.
Intervention periods are therefore, by design, the **later** periods, and any improvement over time
looks like an intervention effect unless the model has **period effects**. Hussey and Hughes (2007):
period as fixed effects, a random intercept per cluster.

| Model | Intervention effect here (truth −2.0) |
|---|---|
| **with period fixed effects** | −1.56 |
| without period effects | **+1.22**: the secular trend, credited to the intervention |

## Few clusters

With 12 wards the Wald interval is too narrow. The final report should use a Kenward-Roger or
Satterthwaite interval (`lmerTest`), or randomization inference over the possible rollout orders. The
Wald interval here is the baseline both languages agree on.

## Three Python details

- The optimizer starts from a moment estimate of the ward variance. From statsmodels' default start,
  L-BFGS reached the zero-variance boundary on CI's Linux runner and the interior optimum on macOS, from
  the same data and package versions. The fit asserts the final ward variance is positive, which is how
  that was caught.

- The intervention's standard error is computed as lme4 reports it: (X′V⁻¹X)⁻¹ conditional on the REML
  variance components. statsmodels' `bse` inverts the joint Hessian of fixed and variance parameters,
  which differs slightly in an unbalanced design (1.8e-5 here).
- statsmodels warns "Random effects covariance is singular" when an optimizer step touches zero; the
  file silences that and asserts the final ward variance is positive.

## The fixture

12 wards, 7 periods, 50 patients per ward-period; two wards switch at each of periods 2 to 7. Secular
trend +1 a period; intervention −2; ward SD 2; noise SD 8.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |

The control file also fits the model with an independent REML written in numpy, and both engines match
it.
