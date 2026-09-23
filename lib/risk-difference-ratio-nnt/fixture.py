"""Seeded fixture for a two-arm trial with a binary outcome.

800 per arm; the outcome occurs in 30% of controls and 20% of the treated. True risk difference
(treated minus control) -0.10, true risk ratio 2/3, true number needed to treat 10.

The committed seed gives a risk-difference interval that excludes zero, so the NNT interval is
finite. When the risk-difference interval crosses zero, the NNT interval runs through infinity and
is reported as two disjoint ranges (Altman 1998); README.md shows how.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261004
PER_ARM = 800
RISK = {1: 0.20, 0: 0.30}
RISK_DIFFERENCE = RISK[1] - RISK[0]
LOG_RISK_RATIO = math.log(RISK[1] / RISK[0])


def rows(seed=SEED, per_arm=PER_ARM):
    rng = random.Random(seed)
    i = 0
    for treated in (1, 0):
        for _ in range(per_arm):
            i += 1
            yield {"id": i, "treated": treated, "outcome": 1 if rng.random() < RISK[treated] else 0}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "treated", "outcome"])
        writer.writeheader()
        writer.writerows(data)
    x1 = sum(r["outcome"] for r in data if r["treated"] == 1)
    x0 = sum(r["outcome"] for r in data if r["treated"] == 0)
    print(f"wrote fixture.csv: treated {x1}/{PER_ARM}, control {x0}/{PER_ARM}")
