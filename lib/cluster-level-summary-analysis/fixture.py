"""Seeded fixture for a cluster-randomized trial with FEW clusters: 24 clinics, 12 per arm.

Each clinic enrols 60 patients. A clinic's infection risk is its arm's risk (0.35 control, 0.20
intervention) plus a clinic-level shift uniform on -0.12 to 0.12, symmetric, so the true difference in
the average clinic risk (intervention minus control) is -0.15. With equal clinic sizes that is also the
patient-level difference.

With 24 clusters, sandwich and GEE standard errors are too small, and an analysis of the 1440
patients as if independent is far too small. Comparing the 24 clinic-level proportions with a t-test on 22
degrees of freedom is simple and valid; that is this entry.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import random

SEED = 20261023
CLUSTERS_PER_ARM = 12
PER_CLUSTER = 60
RISK = {0: 0.35, 1: 0.20}
TRUE_DIFFERENCE = RISK[1] - RISK[0]


def rows(seed=SEED):
    rng = random.Random(seed)
    cluster = 0
    for arm in (1, 0):
        for _ in range(CLUSTERS_PER_ARM):
            cluster += 1
            risk = RISK[arm] + rng.uniform(-0.12, 0.12)
            for _ in range(PER_CLUSTER):
                yield {"cluster": cluster, "arm": arm, "infection": 1 if rng.random() < risk else 0}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["cluster", "arm", "infection"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} patients in {2 * CLUSTERS_PER_ARM} clinics")
