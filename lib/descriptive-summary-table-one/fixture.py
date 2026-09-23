"""Seeded fixture for a descriptive Table 1 in an observational cohort, by exposure group.

Four baseline characteristics with known distributions, so every summary a Table 1 reports is a
property of the generator:

  age       normal, mean 60 unexposed and 65 exposed, SD 11 in both    -> mean (SD), and an SMD
  los_days  log-normal, median exp(1.2), the same in both groups       -> median [Q1, Q3]
  diabetes  binary, 25% unexposed and 35% exposed                      -> n (%), and an SMD
  smoking   never / former / current, 50 / 30 / 20% in both groups     -> n (%) per level

Length of stay is skewed on purpose: its mean sits well above its median, so an implementation
that reported mean (SD) for it, or took quartiles with a different definition, would be caught.

Exposure is assigned with probability 0.5, so the pooled age mean is 62.5 and the pooled SD is
sqrt(11^2 + 2.5^2), because the two group means differ.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20260926
N = 600
AGE_MEAN = {0: 60.0, 1: 65.0}
AGE_SD = 11.0
LOS_MU, LOS_SIGMA = 1.2, 0.6
DIABETES = {0: 0.25, 1: 0.35}
SMOKING = [("never", 0.5), ("former", 0.3), ("current", 0.2)]


def rows(seed=SEED, n=N):
    rng = random.Random(seed)
    for i in range(1, n + 1):
        exposed = 1 if rng.random() < 0.5 else 0
        age = round(rng.gauss(AGE_MEAN[exposed], AGE_SD), 1)
        los = round(math.exp(rng.gauss(LOS_MU, LOS_SIGMA)), 2)
        diabetes = 1 if rng.random() < DIABETES[exposed] else 0
        u, smoking = rng.random(), SMOKING[-1][0]
        cumulative = 0.0
        for level, p in SMOKING:
            cumulative += p
            if u < cumulative:
                smoking = level
                break
        yield {"id": i, "exposed": exposed, "age": age, "los_days": los,
               "diabetes": diabetes, "smoking": smoking}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(data[0]))
        writer.writeheader()
        writer.writerows(data)
    exposed = sum(r["exposed"] for r in data)
    print(f"wrote fixture.csv: {len(data)} rows, {exposed} exposed, {len(data) - exposed} unexposed")
