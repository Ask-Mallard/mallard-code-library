# Risk difference, risk ratio and NNT from a two-arm trial

Absolute and relative effects, reported together. A relative effect alone hides the baseline risk that
decides whether an effect matters to a patient: a risk ratio of 0.67 means an absolute reduction of 10
points at a baseline of 30%, and of 1 point at a baseline of 3%.

| Quantity | Interval |
|---|---|
| Risk difference (treated − control) | Newcombe's hybrid score (method 10), from each arm's Wilson interval |
| Risk ratio | log scale (Katz) |
| Number needed to treat, 1 / |RD| | the risk difference's bounds, inverted |

## The defaults this entry pins

| Language | Default | Pinned |
|---|---|---|
| Python `confint_proportions_2indep` | `method="wald"` for the difference | `method="newcomb"`; `method="log"` for the ratio |
| both | argument order sets the sign | treated first, so the difference is treated minus control |

## When the NNT interval is not an interval

If the risk-difference interval crosses zero, the NNT interval runs through infinity. Report it as two
ranges: benefit from 1/|upper bound| to ∞, and harm from 1/|lower bound| to ∞ (Altman 1998). Both
implementations here **stop** in that case rather than print a single, wrong interval.

## The fixture

800 per arm; risks 0.20 treated and 0.30 control. True difference −0.10, ratio 2/3, NNT 10. Observed:
161/800 against 217/800, difference −0.07, NNT 14.3.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
