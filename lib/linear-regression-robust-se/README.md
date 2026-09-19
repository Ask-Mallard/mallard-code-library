# Linear regression with heteroskedasticity-consistent standard errors

Estimates an adjusted mean difference. The coefficients are ordinary least squares and nothing in
this entry changes them; what it changes is the **standard error**, and therefore every confidence
interval and p value built on it.

## The default this entry exists to pin

Every package has a button labelled "robust" and each one presses a different estimator.

| Language | What the obvious code gives you |
|---|---|
| R | `vcovHC` defaults to **HC3** |
| Python | a bare `.fit()` is the **classical** variance, with no robust default at all |

The two files therefore disagree by default — R's `vcovHC` gives HC3, a bare statsmodels `.fit()`
gives the classical variance — and no output says which one it used. Every file here names **HC3**.

Python is the one worth singling out. R gives you *something* robust when you ask loosely;
statsmodels gives you the classical variance, which on data like this assumes exactly
what the data contradicts. `mod.fit()` and `mod.fit(cov_type="HC3")` differ by a keyword and by
7% of the exposure's standard error.

## What is heteroskedastic here, and why it has to be

The two arms are generated with different spread: sd 0.5 unexposed, 1.7 exposed. That is ordinary
in clinical data, where a treatment that helps some patients a great deal and others not at all
widens the treated arm. Under **constant** variance every robust estimator converges to the
classical one, so a homoskedastic fixture could not tell them apart at all and pinning the option
would be decoration — the same trap the Cox entry's first fixture fell into with tied censoring
times.

Unequal spread does **not** bias the coefficients. It makes the classical interval the wrong width.
A reader who takes the classical interval here is not reading a biased estimate, they are reading
an interval with the wrong width, which is the harder error to notice.

## An intervention that was tried and removed

The covariate was briefly drawn from a heavy-tailed mixture, on the reasoning that HC3 differs from
HC0 only through **leverage** and a tidy normal covariate has almost none. Measured, it moved the
exposure's HC0-to-HC3 gap from 2.31e-4 to 2.14e-4 — that is, not at all, because the exposure is a
binary column whose leverage barely depends on another covariate's tail. It also cost the covariate
estimate 2.7 robust standard errors of accuracy. It is recorded here rather than quietly deleted
because its rationale reads well and is wrong.

## What the agreement check does and does not catch

Measured on this fixture with an independent standard-library least-squares fit:

| Comparison on the exposure's SE | Gap | Against a 1e-4 tolerance |
|---|---|---|
| classical vs HC3 | 6.35e-3 | ~60x — caught with room to spare |
| HC0 vs HC3 | 1.97e-4 | ~2x |
| HC1 vs HC3 | 1.01e-4 | ~1x — no usable margin |

So the cross-engine check catches the mistake that matters — treating the classical SE as robust —
and **cannot be relied on to separate the robust flavours from each other**. That is stated rather
than glossed, because a check credited with more than it does is worse than no check.

The flavour is guarded instead by `must_appear`, which is enforced against the *code* of both
files.

## The recovery tolerance was measured, not guessed

Three hundred alternative seeds of this generator were fitted: the worst of the three misses had a
median of 0.0585, a 95th percentile of 0.1436 and a maximum of 0.2650. The tolerance is 0.25, which
0.3% of seeds would exceed. The committed seed's worst miss is 0.0283.

The two entries before this one set the figure provisionally and were wrong both times.

## What this entry does not do

It is not cluster-robust. A sandwich standard error corrects for unequal variance between
observations that are still **independent**; where observations share a ward, clinic or patient the
correction has to be by cluster, and this entry must not be routed to such a plan. `clustered` is
`false` in `meta.json` for that reason.

It does not check linearity, and it does not make the model right. A correct standard error on a
misspecified mean model is a precise answer to the wrong question.
