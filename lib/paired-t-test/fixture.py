"""Seeded fixture for a paired comparison: the same 60 people measured before and after.

Systolic pressure before is normal around 150; after = before - 5 + noise with SD 6, so the true
mean within-person change (after minus before) is -5 and the true SD of the change is 6.

The before values vary far more between people (SD 15) than the change does, which is exactly when
pairing matters: a two-sample test on the same columns treats the between-person spread as noise and
reports a much wider interval. harness/comparisons_examples.test.py asserts that it does.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import random

SEED = 20260927
N = 60
MEAN_CHANGE = -5.0
SD_CHANGE = 6.0


def rows(seed=SEED, n=N):
    rng = random.Random(seed)
    for i in range(1, n + 1):
        before = rng.gauss(150.0, 15.0)
        after = before + MEAN_CHANGE + rng.gauss(0.0, SD_CHANGE)
        yield {"id": i, "before": round(before, 1), "after": round(after, 1)}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "before", "after"])
        writer.writeheader()
        writer.writerows(data)
    change = sum(r["after"] - r["before"] for r in data) / len(data)
    print(f"wrote fixture.csv: {len(data)} pairs, mean change {change:.3f} (true {MEAN_CHANGE})")
