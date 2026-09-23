"""Seeded fixture for an unordered categorical outcome: discharge destination.

2000 patients, treated with probability 0.5, discharged home, to rehabilitation or to a nursing
facility. Each destination's probability comes from a softmax with HOME as the reference:

  log(P(rehab) / P(home))   = -0.5 + 0.6 x treated
  log(P(nursing) / P(home)) = -1.0 - 0.8 x treated

so the true relative-risk ratios for treatment are exp(0.6) (rehab vs home) and exp(-0.8) (nursing vs
home). The categories have no order, which is why this is multinomial and not ordinal.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261011
N = 2000
REFERENCE = "home"
COEF = {"rehab": (-0.5, 0.6), "nursing": (-1.0, -0.8)}


def rows(seed=SEED, n=N):
    rng = random.Random(seed)
    for i in range(1, n + 1):
        treated = 1 if rng.random() < 0.5 else 0
        weights = {"home": 1.0, **{k: math.exp(a + b * treated) for k, (a, b) in COEF.items()}}
        u, cumulative, total = rng.random(), 0.0, sum(weights.values())
        destination = "nursing"
        for k in ("home", "rehab", "nursing"):
            cumulative += weights[k] / total
            if u < cumulative:
                destination = k
                break
        yield {"id": i, "treated": treated, "destination": destination}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "treated", "destination"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} rows, " +
          ", ".join(f"{k} {sum(r['destination'] == k for r in data)}" for k in ("home", "rehab", "nursing")))
