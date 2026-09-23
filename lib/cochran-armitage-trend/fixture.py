"""Four ordered dose groups (scored 0, 1, 2, 3) of 150 each, with outcome probabilities 0.10, 0.16,
0.22 and 0.28: exactly linear, so the true slope on the probability scale is 0.06 per dose step.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import random

SEED = 20261006
PER_GROUP = 150
RISK = {0: 0.10, 1: 0.16, 2: 0.22, 3: 0.28}
SLOPE = 0.06


def rows(seed=SEED, per_group=PER_GROUP):
    rng = random.Random(seed)
    i = 0
    for dose, p in RISK.items():
        for _ in range(per_group):
            i += 1
            yield {"id": i, "dose": dose, "outcome": 1 if rng.random() < p else 0}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "dose", "outcome"])
        writer.writeheader()
        writer.writerows(data)
    counts = {k: sum(r["outcome"] for r in data if r["dose"] == k) for k in RISK}
    print(f"wrote fixture.csv: events by dose {counts}")
