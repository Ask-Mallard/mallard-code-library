# Instrumental variables by two-stage least squares (Mendelian randomization)

When an exposure is confounded by something unmeasured, a variable that moves the exposure and nothing
else (here a genetic allele score for LDL cholesterol) can recover its effect.

| This fixture | Estimate |
|---|---|
| truth | 0.50 per unit LDL |
| **2SLS** | **0.59** (robust SE 0.076) |
| ordinary regression (confounded) | 1.21 |
| first-stage F | 298 |

## The three assumptions

1. **Relevance**: the instrument moves the exposure. Checkable: the first-stage F should be well above
   10; weak instruments bias 2SLS towards the confounded estimate.
2. **Independence**: nothing causes both the instrument and the outcome.
3. **Exclusion**: the instrument affects the outcome only through the exposure. A pleiotropic variant
   breaks it.

Only the first can be tested from the data. State the others as arguments.

## Do not run the two stages by hand

Regressing the outcome on the fitted exposure gives the right coefficient and the wrong standard error
(28% too large here), because the residuals use the fitted exposure instead of the observed one.
`ivreg` computes it correctly; the Python file writes 2SLS out with the correct residuals. Pinned: HC1.

## Verification

| Engine | Status |
|---|---|
| R (ivreg) | executed in CI |
| Python (written out) | executed in CI |
