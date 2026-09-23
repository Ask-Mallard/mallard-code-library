"""Seeded fixture for a rank comparison of two independent groups with a skewed outcome.

Both groups share one right-skewed (log-normal) shape; the treated group is the same distribution
shifted up by exactly 2 units. Under a pure location shift the Hodges-Lehmann estimator estimates
that shift, so the truth is 2.

A pure shift is also the only case where "the Mann-Whitney test compares medians" is true. When the
shapes differ it compares something else (the probability that a random treated value exceeds a
random control value), which README.md spells out.

40 per group, no ties (values carry six decimals), so R's exact distribution applies.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20260929
PER_GROUP = 40
SHIFT = 2.0


def rows(seed=SEED, per_group=PER_GROUP):
    rng = random.Random(seed)
    i = 0
    for group in (0, 1):
        for _ in range(per_group):
            i += 1
            y = math.exp(rng.gauss(1.0, 0.6)) + SHIFT * group
            yield {"id": i, "group": group, "y": round(y, 6)}


if __name__ == "__main__":
    data = list(rows())
    assert len({r["y"] for r in data}) == len(data), "ties would break the exact distribution"
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "group", "y"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} rows, no ties")
