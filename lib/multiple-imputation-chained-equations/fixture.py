"""Seeded fixture for multiple imputation: a confounder missing at random, depending on the outcome.

1,000 patients. Baseline severity X ~ N(0, 1) confounds treatment and outcome:
logit P(A = 1) = 0.8 X, and Y = 2 + 1.0 A + 1.5 X + noise (SD 1.5).
Severity is missing for about 30% of patients, more often when the outcome is high:
logit P(X missing) = -1.2 + 3.0 (Y - 3) / 2. That is missing at random (it depends on an observed
variable, Y), not completely at random, so dropping the incomplete rows (complete-case analysis) biases
the adjusted effect of A. Multiple imputation with the outcome in the imputation model recovers it.

The adjusted effect of A is 1.0. The fixture writes severity as blank where it is missing.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261115
N = 1000
EFFECT = 1.0


def expit(x):
    return 1.0 / (1.0 + math.exp(-x))


def rows(seed=SEED):
    rng = random.Random(seed)
    for i in range(1, N + 1):
        x = rng.gauss(0.0, 1.0)
        a = 1 if rng.random() < expit(0.8 * x) else 0
        y = 2 + EFFECT * a + 1.5 * x + rng.gauss(0.0, 1.5)
        missing = rng.random() < expit(-1.2 + 3.0 * (y - 3) / 2)
        yield {"id": i, "treated": a, "severity": "" if missing else round(x, 6), "y": round(y, 6)}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "treated", "severity", "y"])
        writer.writeheader()
        writer.writerows(data)
    miss = sum(1 for r in data if r["severity"] == "")
    print(f"wrote fixture.csv: {len(data)} patients, severity missing for {miss}; effect {EFFECT}")
