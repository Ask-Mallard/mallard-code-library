"""Seeded fixture for counts with excess zeros: emergency visits in a year.

800 people, treated with probability 0.5. Thirty percent never visit whatever happens (a structural
zero: they are not at risk, for example because they use another hospital). The rest visit as a Poisson
count with mean 2 (control) or 2 x exp(-0.5) (treated). That is the zero-inflated Poisson model with
true log rate ratio -0.5 in the count part and a structural-zero probability of 0.3, independent of
treatment.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261014
N = 800
ZERO_PROB = 0.3
LOG_MEAN_CONTROL = math.log(2.0)
LOG_RATE_RATIO = -0.5
INFLATION_LOGIT = math.log(ZERO_PROB / (1 - ZERO_PROB))


def rows(seed=SEED, n=N):
    rng = random.Random(seed)
    for i in range(1, n + 1):
        treated = 1 if rng.random() < 0.5 else 0
        structural_zero = rng.random() < ZERO_PROB
        visits = 0
        if not structural_zero:
            rate = math.exp(LOG_MEAN_CONTROL + LOG_RATE_RATIO * treated)
            t = rng.expovariate(rate)
            while t <= 1.0:
                visits += 1
                t += rng.expovariate(rate)
        yield {"id": i, "treated": treated, "visits": visits}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "treated", "visits"])
        writer.writeheader()
        writer.writerows(data)
    zeros = sum(1 for r in data if r["visits"] == 0)
    print(f"wrote fixture.csv: {len(data)} rows, {zeros} zeros ({zeros / len(data):.1%})")
