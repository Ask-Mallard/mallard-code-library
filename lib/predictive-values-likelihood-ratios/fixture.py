"""Seeded fixture for predictive values and likelihood ratios from a diagnostic accuracy study.

1,500 consecutive patients, prevalence 0.2. The test has sensitivity 0.85 and specificity 0.90. So
PPV = 0.85 x 0.2 / (0.85 x 0.2 + 0.10 x 0.8) = 0.68, NPV = 0.90 x 0.8 / (0.90 x 0.8 + 0.15 x 0.2) = 0.96,
LR+ = 0.85 / 0.10 = 8.5 and LR- = 0.15 / 0.90 = 0.1667. At a prevalence of 0.05 the same test has a PPV
of only 0.309: predictive values belong to a population, likelihood ratios to the test.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261121
N = 1500
PREVALENCE = 0.2
SENSITIVITY = 0.85
SPECIFICITY = 0.90


def truth(prevalence=PREVALENCE):
    se, sp, p = SENSITIVITY, SPECIFICITY, prevalence
    return {"ppv": se * p / (se * p + (1 - sp) * (1 - p)), "npv": sp * (1 - p) / (sp * (1 - p) + (1 - se) * p),
            "log_lr_positive": math.log(se / (1 - sp)), "log_lr_negative": math.log((1 - se) / sp)}


def rows(seed=SEED):
    rng = random.Random(seed)
    for i in range(1, N + 1):
        d = 1 if rng.random() < PREVALENCE else 0
        positive = rng.random() < (SENSITIVITY if d else 1 - SPECIFICITY)
        yield {"id": i, "disease": d, "test_positive": 1 if positive else 0}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "disease", "test_positive"])
        writer.writeheader()
        writer.writerows(data)
    t = truth()
    print(f"wrote fixture.csv: {len(data)} patients; PPV {t['ppv']:.4f}, NPV {t['npv']:.4f}, "
          f"LR+ {math.exp(t['log_lr_positive']):.4f}, LR- {math.exp(t['log_lr_negative']):.4f}")
