"""Seeded fixture for a stepped-wedge cluster-randomized trial with a continuous outcome.

12 wards, 7 periods. All wards start in control; at each of periods 2 to 7, two more wards switch to the
intervention, in a random order, and stay switched. Each ward-period measures 50 patients.

The outcome is 50 + 1.0 x (period - 1) + (-2.0) x intervention + ward effect (SD 2) + noise (SD 8). The
secular trend (+1 a period) matters: intervention periods are, by design, the LATER periods, so a model
without period effects credits the trend to the intervention. The true intervention effect is -2.0.

This is the Hussey and Hughes (2007) model: period as fixed effects, a random intercept per ward.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import random

SEED = 20261025
WARDS = 12
PERIODS = 7
PER_CELL = 50
EFFECT = -2.0
TREND = 1.0


def rows(seed=SEED):
    rng = random.Random(seed)
    order = list(range(1, WARDS + 1))
    rng.shuffle(order)
    switch = {w: 2 + i // 2 for i, w in enumerate(order)}   # two wards switch at each of periods 2..7
    ward_effect = {w: rng.gauss(0.0, 2.0) for w in range(1, WARDS + 1)}
    for w in range(1, WARDS + 1):
        for p in range(1, PERIODS + 1):
            treated = 1 if p >= switch[w] else 0
            for _ in range(PER_CELL):
                y = 50 + TREND * (p - 1) + EFFECT * treated + ward_effect[w] + rng.gauss(0.0, 8.0)
                yield {"ward": w, "period": p, "intervention": treated, "y": round(y, 4)}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["ward", "period", "intervention", "y"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} rows, {WARDS} wards x {PERIODS} periods x {PER_CELL}")
