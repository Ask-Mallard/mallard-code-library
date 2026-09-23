"""Seeded fixture for a paired rank comparison: 45 people measured before and after.

The within-person change (after minus before) is drawn from a Laplace distribution centred on +3:
symmetric, but with heavier tails than a normal, so a few changes are large. The signed-rank test
assumes the changes are SYMMETRIC about their centre, and under that assumption the Hodges-Lehmann
pseudo-median estimates the centre, so the truth is 3.

Values carry six decimals so there are no zero changes and no ties, and R's exact distribution
applies. Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import random

SEED = 20260930
N = 45
CENTRE = 3.0
SCALE = 3.0


def rows(seed=SEED, n=N):
    rng = random.Random(seed)
    for i in range(1, n + 1):
        before = rng.gauss(20.0, 5.0)
        change = CENTRE + rng.expovariate(1 / SCALE) - rng.expovariate(1 / SCALE)
        yield {"id": i, "before": round(before, 6), "after": round(before + change, 6)}


if __name__ == "__main__":
    data = list(rows())
    changes = [round(r["after"] - r["before"], 9) for r in data]
    assert 0 not in changes and len({abs(c) for c in changes}) == len(changes), "zeros or ties"
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "before", "after"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} pairs, no zero changes, no tied magnitudes")
