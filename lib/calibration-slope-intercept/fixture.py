"""Seeded fixture for assessing the calibration of a published prediction model in new patients.

2,000 patients in a validation cohort, with two standardized predictors X1, X2 ~ N(0, 1). Their true
risk is logit P(Y = 1) = -1.5 + 0.8 X1 + 0.5 X2. The published model, developed elsewhere and overfitted,
predicts logit = -1.0 + 1.2 X1 + 0.75 X2: its coefficients are 1.5 times too large and its intercept
too high. So the true risk, written in terms of the published linear predictor LP, is
logit P(Y = 1) = -0.8333 + 0.6667 LP: the CALIBRATION SLOPE is 2/3 (predictions too extreme) and the
model over-predicts on average.

`truth()` returns the slope (exactly 2/3), the calibration intercept with the slope fixed at 1
(calibration-in-the-large, found by solving E[expit(a + LP)] = E[P(Y = 1)] by quadrature), and the
expected observed-to-expected ratio. The fixture stores the published model's LP beside the outcome.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261122
N = 2000
TRUE = (-1.5, 0.8, 0.5)
PUBLISHED = (-1.0, 1.2, 0.75)


def expit(x):
    return 1.0 / (1.0 + math.exp(-x))


def rows(seed=SEED):
    rng = random.Random(seed)
    for i in range(1, N + 1):
        x1, x2 = rng.gauss(0.0, 1.0), rng.gauss(0.0, 1.0)
        y = 1 if rng.random() < expit(TRUE[0] + TRUE[1] * x1 + TRUE[2] * x2) else 0
        lp = PUBLISHED[0] + PUBLISHED[1] * x1 + PUBLISHED[2] * x2
        yield {"id": i, "lp": round(lp, 6), "event": y}


def _normal_mean(f, mean, sd, steps=20000):
    lo, hi = mean - 12 * sd, mean + 12 * sd
    h = (hi - lo) / steps
    total = 0.0
    for k in range(steps + 1):
        x = lo + k * h
        w = 1 if k in (0, steps) else (4 if k % 2 else 2)
        total += w * f(x) * math.exp(-0.5 * ((x - mean) / sd) ** 2) / (sd * math.sqrt(2 * math.pi))
    return total * h / 3


def truth():
    # LP is normal: mean -1.0, SD sqrt(1.2^2 + 0.75^2). The true logit is a + b LP with b = 0.8 / 1.2.
    b = TRUE[1] / PUBLISHED[1]
    a = TRUE[0] - b * PUBLISHED[0]
    sd = math.sqrt(PUBLISHED[1] ** 2 + PUBLISHED[2] ** 2)
    observed = _normal_mean(lambda lp: expit(a + b * lp), PUBLISHED[0], sd)
    expected = _normal_mean(expit, PUBLISHED[0], sd)
    lo, hi = -5.0, 5.0
    for _ in range(100):  # bisection for the calibration-in-the-large intercept
        mid = (lo + hi) / 2
        if _normal_mean(lambda lp: expit(mid + lp), PUBLISHED[0], sd) < observed:
            lo = mid
        else:
            hi = mid
    return {"calibration_slope": b, "calibration_in_the_large": (lo + hi) / 2, "observed_expected": observed / expected}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "lp", "event"])
        writer.writeheader()
        writer.writerows(data)
    t = truth()
    print(f"wrote fixture.csv: {len(data)} patients, {sum(r['event'] for r in data)} events; slope "
          f"{t['calibration_slope']:.10f}, CITL {t['calibration_in_the_large']:.10f}, O/E {t['observed_expected']:.10f}")
