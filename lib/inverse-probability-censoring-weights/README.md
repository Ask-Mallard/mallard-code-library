# Inverse probability of censoring weights for informative dropout

When sicker patients drop out more, and more so in one arm, the patients whose outcome is observed do not
represent their arm, and the completers' risk difference is biased, even in a randomized trial. Weighting
each completer by 1 / P(stayed | arm, prognostic factors) lets them stand in for everyone randomized.

| This fixture | Risk difference |
|---|---|
| truth | −0.102 |
| **IPCW** | **−0.083** (SE 0.024) |
| complete-case | −0.046 |

## What it assumes

Missing at random given arm and the modelled factors: no unmeasured cause of both dropout and outcome. The
dropout model here includes the arm-by-factor interaction.

## The standard error

From the stacked estimating equations (the dropout model's score and the two weighted risks), so the
estimation of the weights is included; the weights-as-known SE is reported beside it. With time-to-event
outcomes the same idea gives time-varying weights and a weighted Kaplan-Meier or Cox model (see
`target-trial-clone-censor-weight`).

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
