"""Seeded fixture for a continuous outcome measured repeatedly: four visits per patient.

300 patients, treated with probability 0.5, measured at visits 0, 1, 2 and 3. The outcome is
50 - 1.0 x visit - 1.5 x treated x visit + a patient-level intercept (SD 5) + noise (SD 3), so both arms
start level and the treated decline 1.5 units per visit faster: the true treatment-by-time
interaction is -1.5.

The patient intercept makes a patient's four measurements correlated (intraclass correlation
25 / (25 + 9) = 0.74). Treating the 1200 rows as independent understates the uncertainty of anything
that compares patients; the mixed model's random intercept accounts for it.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import random

SEED = 20261022
N = 300
VISITS = (0, 1, 2, 3)
INTERACTION = -1.5
TIME_SLOPE = -1.0
SD_PATIENT, SD_NOISE = 5.0, 3.0


def rows(seed=SEED, n=N):
    rng = random.Random(seed)
    for i in range(1, n + 1):
        treated = 1 if rng.random() < 0.5 else 0
        u = rng.gauss(0.0, SD_PATIENT)
        for v in VISITS:
            y = 50 + TIME_SLOPE * v + INTERACTION * treated * v + u + rng.gauss(0.0, SD_NOISE)
            yield {"id": i, "treated": treated, "visit": v, "y": round(y, 4)}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "treated", "visit", "y"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} rows for {N} patients")
