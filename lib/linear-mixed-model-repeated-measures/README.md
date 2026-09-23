# Linear mixed model for repeated measurements

Each patient measured at several visits. Their measurements are correlated (here 74% of the variance is
between patients), so a **random intercept** per patient is fitted, and the question "do the arms change
at different rates?" is the **treatment-by-time interaction**.

## The defaults this entry pins

| Choice | Pinned | Why |
|---|---|---|
| REML or ML | REML in both | the default in both, named: they give different variance components and standard errors |
| visit | numeric | a linear trend; as a factor it would estimate a difference at every visit |
| statsmodels optimizer | L-BFGS with `pgtol` and `factr` tightened | at the defaults the patient SD stopped 3.8e-5 short of the exact REML value |

statsmodels warns that `pgtol` and `factr` are "not used by MixedLM.fit". They are passed to the
optimizer, and the fit reaches the exact value, so the Python file silences that one warning on purpose.

## What the random intercept assumes

That every patient changes at the same rate apart from treatment. When individual trajectories differ,
add a random slope (`(1 + visit | id)`); otherwise the interaction's standard error is too small.

## The fixture

300 patients, 4 visits; outcome 50 − 1.0 × visit − 1.5 × treated × visit + patient intercept (SD 5) +
noise (SD 3). True interaction −1.5; observed −1.34.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |

In this balanced design the REML answers have closed forms, and the control file checks both engines
against them.
