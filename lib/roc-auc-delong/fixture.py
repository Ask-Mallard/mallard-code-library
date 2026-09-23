"""Seeded fixture for an ROC analysis: one continuous biomarker against a reference standard.

600 patients, disease prevalence 0.3. The marker is N(0, 1) without disease and N(1.2, 1) with it
(binormal, equal variance), so higher values point to disease. The true area under the ROC curve is
Phi(1.2 / sqrt(2)) = 0.8019, and the threshold that maximizes Youden's index (sensitivity + specificity
- 1) is the midpoint, 0.6.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261119
N = 600
PREVALENCE = 0.3
SHIFT = 1.2


def true_auc():
    return 0.5 * (1 + math.erf(SHIFT / math.sqrt(2) / math.sqrt(2)))


def rows(seed=SEED):
    rng = random.Random(seed)
    for i in range(1, N + 1):
        d = 1 if rng.random() < PREVALENCE else 0
        yield {"id": i, "disease": d, "marker": round(rng.gauss(SHIFT * d, 1.0), 6)}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "disease", "marker"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} patients, {sum(r['disease'] for r in data)} with disease; AUC {true_auc():.10f}")
