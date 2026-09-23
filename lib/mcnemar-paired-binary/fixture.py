"""Seeded fixture for a paired binary outcome: the same 500 people tested before and after.

Positive before with probability 0.30. Positive after with probability 0.80 if positive before and
0.20 if negative before. So 38% are positive after, and the true change in the proportion positive
(after minus before) is 0.08.

Only the DISCORDANT pairs carry information about the change: expected 0.30 x 0.20 = 6% positive ->
negative and 0.70 x 0.20 = 14% negative -> positive. A two-sample chi-square on the same columns
treats before and after as independent groups and ignores that pairing.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import random

SEED = 20261003
N = 500
P_BEFORE = 0.30
P_AFTER_GIVEN = {1: 0.80, 0: 0.20}
TRUE_CHANGE = P_BEFORE * P_AFTER_GIVEN[1] + (1 - P_BEFORE) * P_AFTER_GIVEN[0] - P_BEFORE


def rows(seed=SEED, n=N):
    rng = random.Random(seed)
    for i in range(1, n + 1):
        before = 1 if rng.random() < P_BEFORE else 0
        after = 1 if rng.random() < P_AFTER_GIVEN[before] else 0
        yield {"id": i, "before": before, "after": after}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "before", "after"])
        writer.writeheader()
        writer.writerows(data)
    b = sum(1 for r in data if r["before"] == 1 and r["after"] == 0)
    c = sum(1 for r in data if r["before"] == 0 and r["after"] == 1)
    print(f"wrote fixture.csv: {len(data)} pairs, discordant +/- {b}, -/+ {c}, true change {TRUE_CHANGE:.3f}")
