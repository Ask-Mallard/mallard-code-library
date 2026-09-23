"""Seeded fixture for a cluster trial with few clusters of UNEQUAL size, analysed at the patient level.

12 clinics, 6 per arm, enrolling between 20 and 60 patients each. The outcome is
20 + 4 x treated + 0.1 x (age - 50) + clinic effect (SD 1.5) + noise (SD 5), so the true treatment effect
adjusted for age is 4.

With 12 clusters the ordinary cluster-robust (sandwich) standard error is too small and its normal
reference too optimistic. CR2 (the Bell-McCaffrey bias-reduced sandwich) with Satterthwaite degrees of
freedom is the recommended small-sample correction; unequal cluster sizes are where it differs most
from the ordinary sandwich.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import random

SEED = 20261026
CLINICS_PER_ARM = 6
EFFECT = 4.0
AGE_SLOPE = 0.1


def rows(seed=SEED):
    rng = random.Random(seed)
    clinic = 0
    for treated in (1, 0):
        for _ in range(CLINICS_PER_ARM):
            clinic += 1
            size = rng.randint(20, 60)
            u = rng.gauss(0.0, 1.5)
            for _ in range(size):
                age = round(rng.uniform(30.0, 70.0), 1)
                y = 20 + EFFECT * treated + AGE_SLOPE * (age - 50) + u + rng.gauss(0.0, 5.0)
                yield {"clinic": clinic, "treated": treated, "age": age, "y": round(y, 4)}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["clinic", "treated", "age", "y"])
        writer.writeheader()
        writer.writerows(data)
    sizes = [sum(1 for r in data if r["clinic"] == c) for c in range(1, 2 * CLINICS_PER_ARM + 1)]
    print(f"wrote fixture.csv: {len(data)} patients in {len(sizes)} clinics, sizes {min(sizes)} to {max(sizes)}")
