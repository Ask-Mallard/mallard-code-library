"""Seeded fixture for instrumental variables: a Mendelian-randomization-style analysis.

An allele score Z (the count of risk alleles over four variants, each Binomial(2, 0.3)) raises LDL
cholesterol X. An unmeasured factor U raises both LDL and the outcome Y, so the ordinary regression of
Y on X is confounded:

X = 3 + 0.25 Z + 0.8 U + noise (SD 0.6)
Y = 1 + 0.5 X + 1.0 U + noise (SD 1)

The causal effect of X on Y is 0.5. Z affects Y only through X and is independent of U, which is what
makes it an instrument; the three IV assumptions hold by construction here and cannot be checked from
the data in general.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import random

SEED = 20261109
N = 3000
EFFECT = 0.5


def rows(seed=SEED):
    rng = random.Random(seed)
    for i in range(1, N + 1):
        z = sum(1 for _ in range(8) if rng.random() < 0.3)  # four variants, two alleles each
        u = rng.gauss(0.0, 1.0)
        x = 3 + 0.25 * z + 0.8 * u + rng.gauss(0.0, 0.6)
        y = 1 + EFFECT * x + 1.0 * u + rng.gauss(0.0, 1.0)
        yield {"id": i, "allele_score": z, "ldl": round(x, 6), "y": round(y, 6)}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "allele_score", "ldl", "y"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} participants; causal effect {EFFECT}")
