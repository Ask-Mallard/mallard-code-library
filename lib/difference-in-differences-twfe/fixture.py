"""Seeded fixture for difference-in-differences: a policy adopted by some hospitals at one time.

40 hospitals observed for 8 quarters. 20 adopt a discharge policy at quarter 5 (all at once: not
staggered); 20 never do. The outcome is mean length of stay per hospital-quarter:
y = 6 + 3 [adopter] + hospital effect (SD 1.5) + 0.4 x quarter + effect x [adopter and quarter >= 5]
    + noise (SD 0.8), with effect = -1.2.

Adopters have longer stays to begin with (the +3), and stays drift upwards over time for everyone (the
0.4 a quarter), so neither the before-after change in adopters nor the after-period comparison between
groups estimates the effect. Parallel trends hold by construction.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import random

SEED = 20261106
HOSPITALS = 40
QUARTERS = 8
ADOPT_AT = 5
EFFECT = -1.2


def rows(seed=SEED):
    rng = random.Random(seed)
    for h in range(1, HOSPITALS + 1):
        adopter = 1 if h <= HOSPITALS // 2 else 0
        u = rng.gauss(0.0, 1.5)
        for q in range(1, QUARTERS + 1):
            post = 1 if q >= ADOPT_AT else 0
            y = 6 + 3 * adopter + u + 0.4 * q + EFFECT * adopter * post + rng.gauss(0.0, 0.8)
            yield {"hospital": h, "quarter": q, "adopter": adopter, "post": post,
                   "policy": adopter * post, "los": round(y, 6)}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["hospital", "quarter", "adopter", "post", "policy", "los"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {HOSPITALS} hospitals x {QUARTERS} quarters; effect {EFFECT}")
