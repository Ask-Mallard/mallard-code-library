"""Seeded fixture for an ordered outcome: none, mild, moderate, severe.

500 people, treated with probability 0.5. A latent severity = 0.8 x treated + logistic noise is cut at
-1, 0.5 and 2 into four ordered categories. That is exactly the proportional-odds (cumulative logit)
model, with true log odds ratio 0.8 for being in a HIGHER category and thresholds -1, 0.5, 2.

The category names are words whose ALPHABETICAL order (mild, moderate, none, severe) is not their
clinical order, on purpose: a model handed the column as unordered text sorts it alphabetically and
fits a different, meaningless ordering without any warning.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261010
N = 500
LOG_OR = 0.8
CUTS = (-1.0, 0.5, 2.0)
LEVELS = ("none", "mild", "moderate", "severe")


def rows(seed=SEED, n=N):
    rng = random.Random(seed)
    for i in range(1, n + 1):
        treated = 1 if rng.random() < 0.5 else 0
        u = rng.random()
        latent = LOG_OR * treated + math.log(u / (1 - u))
        level = LEVELS[sum(latent > c for c in CUTS)]
        yield {"id": i, "treated": treated, "severity": level}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "treated", "severity"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} rows, " +
          ", ".join(f"{lv} {sum(r['severity'] == lv for r in data)}" for lv in LEVELS))
