"""Seeded fixture for decision curve analysis: the net benefit of acting on a prediction model.

2,000 patients. A well-calibrated risk model gives each patient their true risk, logit p = -1.5 + 1.2 X
with X ~ N(0, 1); the outcome is Bernoulli(p). A clinician treats when the predicted risk exceeds a
threshold t, where t expresses how many false positives are worth one true positive (odds t / (1 - t)).

Net benefit at threshold t = TP/n - FP/n x t / (1 - t). `truth()` integrates it for the model and for
"treat everyone" at t = 0.1, 0.2 and 0.3 from the data-generating process.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261126
N = 2000
INTERCEPT, SLOPE = -1.5, 1.2
THRESHOLDS = (0.1, 0.2, 0.3)


def expit(x):
    return 1.0 / (1.0 + math.exp(-x))


def rows(seed=SEED):
    rng = random.Random(seed)
    for i in range(1, N + 1):
        p = expit(INTERCEPT + SLOPE * rng.gauss(0.0, 1.0))
        yield {"id": i, "predicted_risk": round(p, 6), "event": 1 if rng.random() < p else 0}


def _normal_mean(f, steps=40000):
    lo, hi = -12.0, 12.0
    h = (hi - lo) / steps
    total = 0.0
    for k in range(steps + 1):
        x = lo + k * h
        w = 1 if k in (0, steps) else (4 if k % 2 else 2)
        total += w * f(x) * math.exp(-0.5 * x * x) / math.sqrt(2 * math.pi)
    return total * h / 3


def truth():
    out = {}
    prevalence = _normal_mean(lambda x: expit(INTERCEPT + SLOPE * x))
    for t in THRESHOLDS:
        # The model treats when p > t, i.e. when x > (logit t - intercept) / slope; integrate on that side.
        cut = (math.log(t / (1 - t)) - INTERCEPT) / SLOPE
        above = lambda f: _normal_mean(lambda x: f(x) if x > cut else 0.0)
        tp = above(lambda x: expit(INTERCEPT + SLOPE * x))
        fp = above(lambda x: 1 - expit(INTERCEPT + SLOPE * x))
        key = f"{t:.1f}".replace(".", "")
        out[f"net_benefit_model_t{key}"] = tp - fp * t / (1 - t)
        out[f"net_benefit_all_t{key}"] = prevalence - (1 - prevalence) * t / (1 - t)
    return out


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "predicted_risk", "event"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} patients; " + ", ".join(f"{k} {v:.6f}" for k, v in truth().items()))
