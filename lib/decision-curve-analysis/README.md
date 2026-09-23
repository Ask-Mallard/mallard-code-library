# Decision curve analysis: net benefit of a prediction model

A model can discriminate well and still not improve decisions. Decision curve analysis (Vickers and
Elkin 2006) asks, at a threshold risk t above which a clinician would treat:

net benefit = TP/n − FP/n × t / (1 − t)

and compares it with **treat everyone** (prevalence − (1 − prevalence) × t / (1 − t)) and **treat no one**
(0). The model is useful at a threshold only if it beats both.

| This fixture | Model: truth / estimate | Treat all: truth / estimate |
|---|---|---|
| t = 0.1 | 0.163 / 0.154 | 0.149 / 0.143 |
| t = 0.2 | 0.111 / 0.094 | 0.043 / 0.036 |
| t = 0.3 | 0.073 / 0.068 | −0.094 / −0.102 |

Choose the range of thresholds from the clinical decision beforehand and plot the whole curve over it.
Net benefit uses the predicted risks as given, so a miscalibrated model is penalized, as intended.
Intervals, if wanted, come from a bootstrap.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
