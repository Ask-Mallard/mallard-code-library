"""Seeded fixture for a single proportion: one binary outcome in a simple random sample.

The true proportion is 0.12, well away from 0.5 on purpose. Near 0.5 the Wald, Wilson and
Clopper-Pearson intervals nearly coincide, so an implementation that silently fell back to Wald
would pass. At 0.12 with n = 400 they separate by more than the agreement tolerance, which
harness/descriptive_examples.test.py asserts.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV, and the
controls regenerate it byte for byte.
"""

import csv
import random

SEED = 20260923
N = 400
PROPORTION = 0.12


def rows(seed=SEED, n=N, p=PROPORTION):
    rng = random.Random(seed)
    for i in range(1, n + 1):
        yield {"id": i, "outcome": 1 if rng.random() < p else 0}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "outcome"])
        writer.writeheader()
        writer.writerows(data)
    events = sum(r["outcome"] for r in data)
    print(f"wrote fixture.csv: {len(data)} rows, {events} events, "
          f"observed {events / len(data):.4f} (true {PROPORTION})")
