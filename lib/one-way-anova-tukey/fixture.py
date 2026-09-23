"""Seeded fixture for comparing a continuous outcome across three groups.

Three arms of 100, coded 1, 2 and 3, with true means 10, 12 and 15 and a common SD of 2.5. The pairwise
mean differences are therefore 2 (arm 2 - arm 1), 5 (3 - 1) and 3 (3 - 2), and the pooled SD is 2.5.

The arms are coded as NUMBERS on purpose. A model that reads the arm column as numeric fits one
straight-line slope (1 degree of freedom) instead of three group means (2 degrees of freedom), with
no warning. The unequal spacing of the means (2 then 3) makes that fit visibly different.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import random

SEED = 20260928
PER_ARM = 100
MEANS = {1: 10.0, 2: 12.0, 3: 15.0}
SD = 2.5


def rows(seed=SEED, per_arm=PER_ARM):
    rng = random.Random(seed)
    i = 0
    for arm, mean in MEANS.items():
        for _ in range(per_arm):
            i += 1
            yield {"id": i, "arm": arm, "y": round(rng.gauss(mean, SD), 2)}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "arm", "y"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} rows in {len(MEANS)} arms")
