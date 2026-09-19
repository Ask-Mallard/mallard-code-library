# Modified Poisson regression for an adjusted risk ratio

Logistic regression gives an **odds ratio**. For a common outcome that is not a risk ratio, and it
is not collapsible, so it changes when covariates are added even with no confounding. When the
question asks about risk, fit the risk.

Modified Poisson (Zou 2004) fits a log-link Poisson model to binary data and reads the coefficient
as a log risk ratio. The Poisson variance assumption is wrong by construction, which is why a
**robust variance is mandatory rather than optional**.

## The default this entry exists to pin

`sandwich::vcovHC` in R defaults to **HC3**. Zou's method is the ordinary sandwich, **HC0**. A file
that calls `vcovHC(fit)` without `type` reports a different standard error from the method it cites,
and nothing in the output says so.

statsmodels has the mirror-image trap: GLM defaults to the **model-based** covariance, so
`.fit()` without `cov_type` returns intervals from a variance assumption that does not hold.

Both files pin HC0, and the agreement check on the standard error is what keeps them pinned. This
is the failure the library was built for: two correct-looking implementations, different answers,
nothing visible in either file alone.

## The fixture

`fixture.py` writes 3000 rows with a 12.3% event rate. The outcome is drawn from a log-link risk
model, so the true risk ratio is exactly 1.5 and is a property of the data, not of the analysis.

The generator refuses to run if any fitted probability reaches 1, since a log link does not
constrain it and data no risk model can represent would make the entry meaningless.

Sanity check needing no model: the crude risk ratio in the generated data is **1.4975** against a
true 1.5000.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
