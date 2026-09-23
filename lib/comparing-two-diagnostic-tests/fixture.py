"""Seeded fixture for comparing two diagnostic tests measured on the same patients.

2,000 patients, prevalence 0.3. Two markers, A and B, are measured on everyone. Within each disease
group they are bivariate normal with correlation 0.5 and unit variances; with disease, A is shifted by
1.2 and B by 0.7. True AUCs: Phi(1.2 / sqrt 2) = 0.8019 for A and Phi(0.7 / sqrt 2) = 0.6897 for B, a
difference of 0.1122. At a prespecified common threshold of 0.6 for both markers,
the sensitivities are Phi(0.6) = 0.7257 and Phi(0.1) = 0.5398, a difference of 0.1859.

Because the same patients have both tests, the AUCs and the sensitivities are correlated, and the
comparison must be paired.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261120
N = 2000
PREVALENCE = 0.3
SHIFT_A, SHIFT_B = 1.2, 0.7
RHO = 0.5
THRESHOLD_A, THRESHOLD_B = 0.6, 0.6


def phi(x):
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


def truth():
    auc_a, auc_b = phi(SHIFT_A / math.sqrt(2)), phi(SHIFT_B / math.sqrt(2))
    sens_a, sens_b = phi(SHIFT_A - THRESHOLD_A), phi(SHIFT_B - THRESHOLD_B)
    return {"auc_a": auc_a, "auc_b": auc_b, "auc_difference": auc_a - auc_b,
            "sensitivity_difference": sens_a - sens_b}


def rows(seed=SEED):
    rng = random.Random(seed)
    for i in range(1, N + 1):
        d = 1 if rng.random() < PREVALENCE else 0
        z1, z2 = rng.gauss(0.0, 1.0), rng.gauss(0.0, 1.0)
        a = SHIFT_A * d + z1
        b = SHIFT_B * d + RHO * z1 + math.sqrt(1 - RHO ** 2) * z2
        yield {"id": i, "disease": d, "marker_a": round(a, 6), "marker_b": round(b, 6)}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "disease", "marker_a", "marker_b"])
        writer.writeheader()
        writer.writerows(data)
    t = truth()
    print(f"wrote fixture.csv: {len(data)} patients; AUC difference {t['auc_difference']:.10f}, "
          f"sensitivity difference {t['sensitivity_difference']:.10f}")
