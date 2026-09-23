"""Seeded fixture for a prevalence estimated from a cluster sample.

40 clinics are sampled with equal probability from 400, and every eligible patient in a sampled
clinic is included. Patients in one clinic resemble each other, so the rows are not independent.

CLINIC SIZE AND PREVALENCE ARE RELATED ON PURPOSE. Small clinics (10 patients) run at about 30%,
medium (30) at 20%, large (60) at 10%. Two different numbers can then be called "the prevalence":

  * the patient-level prevalence, what share of PATIENTS have the condition:
      (10 * 0.30 + 30 * 0.20 + 60 * 0.10) / (10 + 30 + 60) = 0.15
  * the average of the clinic prevalences, what share a typical CLINIC has:
      (0.30 + 0.20 + 0.10) / 3 = 0.20

This entry estimates the first, which is almost always the question asked. Averaging the clinic
percentages lands near 0.20 and fails recovery, which harness/descriptive_examples.test.py asserts.
When sizes and prevalence are unrelated, both numbers coincide and the mistake cannot be seen.

Each clinic's own prevalence is its size class's value plus uniform noise of +/- 0.05, symmetric,
so the patient-level truth stays 0.15. Stdlib only.
"""

import csv
import random

SEED = 20260925
CLUSTERS = 40
POPULATION_CLUSTERS = 400
SIZE_CLASSES = [(10, 0.30), (30, 0.20), (60, 0.10)]
TRUE_PREVALENCE = sum(m * p for m, p in SIZE_CLASSES) / sum(m for m, _ in SIZE_CLASSES)


def rows(seed=SEED, clusters=CLUSTERS):
    rng = random.Random(seed)
    weight = POPULATION_CLUSTERS / clusters
    for c in range(1, clusters + 1):
        size, base = SIZE_CLASSES[rng.randrange(len(SIZE_CLASSES))]
        risk = base + rng.uniform(-0.05, 0.05)
        for _ in range(size):
            yield {"cluster": c, "weight": weight, "outcome": 1 if rng.random() < risk else 0}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["cluster", "weight", "outcome"])
        writer.writeheader()
        writer.writerows(data)
    n = len(data)
    events = sum(r["outcome"] for r in data)
    print(f"wrote fixture.csv: {CLUSTERS} clusters, {n} patients, "
          f"patient-level prevalence {events / n:.4f} (true {TRUE_PREVALENCE:.4f})")
