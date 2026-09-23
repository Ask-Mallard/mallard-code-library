"""Seeded fixture for an adjusted risk difference from a binary outcome.

1000 people, treated with probability 0.5, age uniform from 40 to 80. The outcome risk is LINEAR on
the probability scale: 0.10 + 0.08 x treated + 0.004 x (age - 60). So the true adjusted risk
difference for treatment is 0.08, and each year of age adds 0.004. Risks run from 0.02 to 0.26, inside
(0, 1), so an identity-link binomial model is correctly specified and can converge.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import random

SEED = 20261008
N = 1000
BASE, RISK_DIFFERENCE, AGE_SLOPE = 0.10, 0.08, 0.004


def rows(seed=SEED, n=N):
    rng = random.Random(seed)
    for i in range(1, n + 1):
        treated = 1 if rng.random() < 0.5 else 0
        age = round(rng.uniform(40.0, 80.0), 1)
        risk = BASE + RISK_DIFFERENCE * treated + AGE_SLOPE * (age - 60)
        yield {"id": i, "treated": treated, "age": age, "outcome": 1 if rng.random() < risk else 0}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "treated", "age", "outcome"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} rows, {sum(r['outcome'] for r in data)} events")
