# Causal mediation with a binary outcome

How much of a treatment's effect on a binary outcome runs through a mediator? With a binary outcome,
**Baron-Kenny and product-of-coefficients are invalid**: the odds ratio is non-collapsible, so the change
in the treatment coefficient when the mediator is added is not an indirect effect, and neither method
has a term for a treatment-by-mediator interaction.

## The counterfactual definitions

With E[Y(a, M(a*))] the risk under treatment a with the mediator as it would be under a*:

- natural direct effect: E[Y(1, M(0))] − E[Y(0, M(0))]
- natural indirect effect: E[Y(1, M(1))] − E[Y(1, M(0))]

They sum to the total effect. The **mediation formula** computes each from a mediator model and an
outcome model with the interaction, averaged over the covariates (VanderWeele 2015). With a binary
mediator the sum over its two values is exact, so the result is reproducible. R's `mediation` package
computes the same quantities by quasi-Bayesian simulation.

| This fixture (risk difference) | Truth | Estimate |
|---|---|---|
| natural direct effect | 0.102 | 0.080 (SE 0.015) |
| natural indirect effect | 0.079 | 0.086 (SE 0.008) |

On the log-odds scale, over 100 seeds, the counterfactual indirect effect averages 0.324 against a truth
of 0.323; the difference method averages 0.257.

## Assumptions

No unmeasured confounding of treatment-outcome, mediator-outcome and treatment-mediator, and no
mediator-outcome confounder that treatment affects. Plan a sensitivity analysis for the mediator-outcome
assumption. **Power for mediation needs its own simulation.**

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |

Both write the formula and the stacked sandwich out; the control file rebuilds the sandwich numerically.
