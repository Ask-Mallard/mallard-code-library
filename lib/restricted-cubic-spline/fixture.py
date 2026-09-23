"""Seeded fixture for a nonlinear relationship: systolic pressure across age.

600 adults, age uniform 20 to 90. The true mean pressure is itself a restricted (natural) cubic
spline in age, with knots at 30, 45, 60, 70 and 80: cubic between the outer knots, linear beyond them.
Because the truth lies exactly in the fitted model's function space, every prediction and contrast has
an exact true value (`curve` below), and recovery checks the fit, not an approximation error.

Noise is normal with SD 10. Stdlib only. `rows(seed)` is the generator; running the file writes the
committed CSV. The coefficients were chosen by fitting this spline to a sigmoid rising from about 115
to 150 mmHg around age 58, then rounded; the rounded spline IS the truth.
"""

import csv
import random

SEED = 20261013
N = 600
KNOTS = (30.0, 45.0, 60.0, 70.0, 80.0)
INTERCEPT, LINEAR = 113.508, 0.07441
SPLINE = (0.0008328, -0.001684, -0.0007937)   # on Harrell's three nonlinear terms
SD = 10.0


def harrell_terms(x, t=KNOTS):
    """Harrell's restricted cubic spline terms (unnormalized): linear beyond the outer knots."""
    k = len(t)
    pos = lambda v: max(v, 0.0) ** 3
    return [pos(x - t[j]) - pos(x - t[k - 2]) * (t[k - 1] - t[j]) / (t[k - 1] - t[k - 2])
            + pos(x - t[k - 1]) * (t[k - 2] - t[j]) / (t[k - 1] - t[k - 2]) for j in range(k - 2)]


def curve(age):
    return INTERCEPT + LINEAR * age + sum(c * h for c, h in zip(SPLINE, harrell_terms(age)))


def rows(seed=SEED, n=N):
    rng = random.Random(seed)
    for i in range(1, n + 1):
        age = round(rng.uniform(20.0, 90.0), 2)
        yield {"id": i, "age": age, "sbp": round(curve(age) + rng.gauss(0.0, SD), 3)}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "age", "sbp"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} rows; true curve at 40/50/70: "
          f"{curve(40):.2f} / {curve(50):.2f} / {curve(70):.2f}")
