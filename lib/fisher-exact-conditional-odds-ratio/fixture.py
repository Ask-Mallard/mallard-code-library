"""Seeded fixture for a two-group comparison of a binary outcome, analysed with Fisher's exact test.

Two groups of 120; the outcome occurs with probability 0.5 in the exposed and 0.2 in the unexposed,
so the true odds ratio is (0.5 / 0.5) / (0.2 / 0.8) = 4, a log odds ratio of log(4).

Fisher's exact test is valid at any size and is the usual choice when expected cell counts are small;
120 per group keeps the recovery check able to tell the odds ratio of 4 from 1. At 30 and 60 per
group the 99th percentile misses (1.34 and 1.09 on the log scale) reached log(4) itself.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261002
PER_GROUP = 120
RISK = {1: 0.5, 0: 0.2}
LOG_OR = math.log((RISK[1] / (1 - RISK[1])) / (RISK[0] / (1 - RISK[0])))


def rows(seed=SEED, per_group=PER_GROUP):
    rng = random.Random(seed)
    i = 0
    for exposed in (1, 0):
        for _ in range(per_group):
            i += 1
            yield {"id": i, "exposed": exposed, "outcome": 1 if rng.random() < RISK[exposed] else 0}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "exposed", "outcome"])
        writer.writeheader()
        writer.writerows(data)
    cells = {(e, o): sum(1 for r in data if r["exposed"] == e and r["outcome"] == o) for e in (1, 0) for o in (1, 0)}
    print(f"wrote fixture.csv: {len(data)} rows, 2x2 {cells}")
