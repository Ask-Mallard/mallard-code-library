"""Seeded fixture for median regression of a skewed outcome, length of stay.

401 patients, treated with probability 0.5, age uniform 40 to 80. The count is odd so the median
linear program has a unique solution (with an even count, a median is an interval and two exact
solvers may return different points of it). Length of stay =
3 + 1.5 x treated + 0.05 x (age - 60) + s x (exp(N(0, 0.8)) - 1), where the error has MEDIAN zero
(exp of a zero-median normal has median 1) and scale s = 1 for controls and 1.8 for the treated.

So the conditional MEDIAN is exactly linear: the true median effect of treatment is 1.5 and of each
year of age 0.05. The conditional MEAN is not the same: the error's mean is exp(0.32) - 1 = 0.377,
scaled by 1.8 in the treated, so the mean effect of treatment is 1.5 + 0.8 x 0.377 = 1.80. Least
squares estimates that mean effect; median regression estimates 1.5. The spread also differs by arm,
which is why the standard errors are the local ('nid') kind rather than ones that assume one error
distribution for everybody.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261012
N = 401
MEDIAN_EFFECT, AGE_SLOPE = 1.5, 0.05
SCALE = {0: 1.0, 1: 1.8}
MEAN_EFFECT = MEDIAN_EFFECT + (SCALE[1] - SCALE[0]) * (math.exp(0.8 ** 2 / 2) - 1)


def rows(seed=SEED, n=N):
    rng = random.Random(seed)
    for i in range(1, n + 1):
        treated = 1 if rng.random() < 0.5 else 0
        age = round(rng.uniform(40.0, 80.0), 1)
        error = SCALE[treated] * (math.exp(rng.gauss(0.0, 0.8)) - 1)
        los = 3 + MEDIAN_EFFECT * treated + AGE_SLOPE * (age - 60) + error
        yield {"id": i, "treated": treated, "age": age, "los_days": round(los, 6)}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "treated", "age", "los_days"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} rows (true median effect {MEDIAN_EFFECT}, mean effect {MEAN_EFFECT:.3f})")
