# Proportional-odds regression for an ordered outcome

For an outcome with ordered categories (none, mild, moderate, severe; a modified Rankin score). One
odds ratio summarizes the shift toward **higher** categories across every cut point together.

## The assumption

The same odds ratio at every cut: none vs mild-or-worse, up to severe vs the rest. Check it by fitting
the binary splits separately, or a partial proportional-odds model, and state it.

## The defaults this entry pins

| Trap | Consequence | Pinned |
|---|---|---|
| category order taken from the text | alphabetical: mild, moderate, **none**, severe; the model fits that ordering silently (log odds ratio 0.35 here instead of 0.98) | the levels written out in both files |
| direction convention | R `polr` and statsmodels model higher categories as positive; SAS `PROC LOGISTIC` by default models lower ones, inverting the odds ratio | stated in both files |
| convergence | statsmodels' BFGS did not converge at its defaults; the two differed by 6.5e-5 | Newton in Python, `reltol = 1e-14` in R |

## The fixture

500 people; latent severity = 0.8 × treated + logistic noise, cut at −1, 0.5 and 2. That is exactly
the proportional-odds model: true log odds ratio 0.8. Observed 0.98.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
