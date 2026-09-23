"""Seeded fixture for an interrupted time series: monthly rates before and after an intervention.

72 months; the intervention starts at month 37 (36 before, 36 after). The monthly outcome follows
y = 50 + 0.2 t - 4 [t >= 37] - 0.3 (t - 36)+ + e_t, with AR(1) errors e_t = 0.5 e_(t-1) + N(0, 0.6^2):
a level drop of 4 and a slope change of -0.3 a month. Neighbouring months are correlated, so ordinary
least squares standard errors are too small.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261107
MONTHS = 72
START = 37
LEVEL_CHANGE = -4.0
SLOPE_CHANGE = -0.3
RHO = 0.5
SD = 0.6


def rows(seed=SEED):
    rng = random.Random(seed)
    e = rng.gauss(0.0, SD / math.sqrt(1 - RHO ** 2))  # stationary start
    for t in range(1, MONTHS + 1):
        if t > 1:
            e = RHO * e + rng.gauss(0.0, SD)
        post = 1 if t >= START else 0
        since = max(t - (START - 1), 0)
        y = 50 + 0.2 * t + LEVEL_CHANGE * post + SLOPE_CHANGE * since + e
        yield {"month": t, "post": post, "months_since": since, "rate": round(y, 6)}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["month", "post", "months_since", "rate"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {MONTHS} months; level {LEVEL_CHANGE}, slope {SLOPE_CHANGE}, AR(1) {RHO}")
