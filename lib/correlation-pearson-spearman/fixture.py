"""Seeded fixture for the correlation between two continuous measurements.

BMI and systolic pressure drawn as a bivariate normal pair with correlation 0.5 (standard normals
rescaled to BMI 27 +/- 4 and SBP 130 +/- 15; rescaling does not change a correlation). For a bivariate
normal, Spearman's rank correlation has a known value too: (6 / pi) * asin(rho / 2) = 0.4826.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261005
N = 150
RHO = 0.5
SPEARMAN = 6 / math.pi * math.asin(RHO / 2)


def rows(seed=SEED, n=N):
    rng = random.Random(seed)
    for i in range(1, n + 1):
        u = rng.gauss(0.0, 1.0)
        v = RHO * u + math.sqrt(1 - RHO ** 2) * rng.gauss(0.0, 1.0)
        yield {"id": i, "bmi": round(27 + 4 * u, 6), "sbp": round(130 + 15 * v, 6)}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "bmi", "sbp"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} rows (true Pearson {RHO}, true Spearman {SPEARMAN:.4f})")
