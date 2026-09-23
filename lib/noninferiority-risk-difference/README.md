# Noninferiority analysis of a risk difference against a margin

A noninferiority trial asks whether a new treatment is **not worse than the standard by more than a
prespecified margin**. The conclusion comes from the confidence interval: noninferior if the lower limit
of the 95% interval for (new − standard) lies above the margin. "No significant difference" is not
noninferiority.

| This fixture (margin −0.10) | Difference | 95% CI (Newcombe) | Noninferior? |
|---|---|---|---|
| intention-to-treat | −0.038 | −0.091 to 0.017 | yes |
| per-protocol | −0.050 | −0.106 to 0.006 | **no** |

## ITT is not conservative here

In a superiority trial intention-to-treat is conservative. In a noninferiority trial it is not:
patients who switch treatments blur the arms together and pull the difference towards 0, towards a
noninferiority conclusion. Report both populations; the conclusion should hold in both. On this fixture
it does not.

## The defaults this entry pins

- The **margin** is set in the protocol from clinical and historical grounds, never after the data.
- **Newcombe's hybrid score interval** (method 10) for the difference; the Wald limit is reported beside
  it.

## Verification

| Engine | Status |
|---|---|
| R (Newcombe written out) | executed in CI |
| Python (statsmodels `confint_proportions_2indep`) | executed in CI |
