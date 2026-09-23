"""Seeded fixture for a sharp regression discontinuity: treatment assigned by a threshold on a score.

Patients with a risk score at or above a cut-off (centred at 0) receive an intervention; those below do
not. The outcome rises smoothly and non-linearly with the score, and jumps by 1.5 at the cut-off:
y = 5 + 2 x + 2 x^3 + 1.5 [x >= 0] + noise (SD 1), with x ~ Uniform(-1, 1).

The effect identified is the jump AT the cut-off, for patients near it. The curvature makes a global
straight line on each side misestimate the jump, which is why local linear fits near the cut-off are
the standard.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import random

SEED = 20261105
N = 2000
JUMP = 1.5


def rows(seed=SEED):
    rng = random.Random(seed)
    for i in range(1, N + 1):
        x = rng.uniform(-1.0, 1.0)
        treated = 1 if x >= 0 else 0
        y = 5 + 2 * x + 2 * x ** 3 + JUMP * treated + rng.gauss(0.0, 1.0)
        yield {"id": i, "score": round(x, 6), "treated": treated, "y": round(y, 6)}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "score", "treated", "y"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} patients, {sum(r['treated'] for r in data)} above the cut-off; jump {JUMP}")
