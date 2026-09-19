# Conditional logistic regression for a matched case-control study

Estimates the exposure odds ratio **conditional on the matched sets**. The matching factor is
eliminated from the likelihood rather than adjusted for, which is why it must not also appear as a
covariate: adding it would condition on the thing already conditioned on.

## What this does not do

**It does not size the study.** A matched or nested case-control has no closed form in Mallard's
sizing catalog, and this entry does not supply one. Borrowing the nearest available formula is a
documented error: `paired` is a continuous mean-difference formula that takes a delta and the SD of
differences, and using it here reports pairs of a quantity nobody measured. The honest route is the
discordant-set method (Dupont 1988) named in prose, with `epiR::epi.sscc` as the reference for
anyone with R who wants to implement it properly.

## The fixture

`fixture.py` writes `fixture.csv`: 400 matched sets, 1 case and 2 controls each, 1200 rows.

The case within each set is drawn with probability proportional to `exp(b1*exposed + b2*covariate)`,
which **is** the conditional likelihood that `clogit` maximises. So the true betas are the fitted
model's own parameters rather than those of some process that merely resembles it, and recovering
them means something.

A large stratum-level term (SD 2.0) is generated and never used in the likelihood. That is
deliberate: matching removes it, so it must not move the estimate. An implementation that silently
fits an *unconditional* model is biased by it and fails the recovery check.

Sanity check available without running any model: the crude odds ratio in the generated data is
**1.983** against a true 2.0.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |

Two separate claims, per `expected.json`: **agreement** (the executed engines produce the same
estimates on identical rows, tolerance 1e-6) and **recovery** (those estimates are near the
parameters the fixture used). Agreement alone would pass two identical mistakes; recovery alone
would pass a default mismatch smaller than sampling error.
