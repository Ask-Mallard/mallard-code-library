"""Seeded fixture for comparing a skewed outcome across three groups.

Three arms of 80, each log-normal with sigma 0.5 and log-scale means 1.0, 1.3 and 1.6, so the true
group medians are exp(1.0), exp(1.3) and exp(1.6). The groups share one shape and differ in
location, which is the situation where the Kruskal-Wallis test reads as a comparison of medians.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261001
PER_ARM = 80
LOG_MEANS = {1: 1.0, 2: 1.3, 3: 1.6}
SIGMA = 0.5


def rows(seed=SEED, per_arm=PER_ARM):
    rng = random.Random(seed)
    i = 0
    for arm, mu in LOG_MEANS.items():
        for _ in range(per_arm):
            i += 1
            yield {"id": i, "arm": arm, "y": round(math.exp(rng.gauss(mu, SIGMA)), 6)}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "arm", "y"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} rows in {len(LOG_MEANS)} arms")
