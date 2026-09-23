"""Seeded fixture for a two-arm randomized trial with a continuous outcome measured at baseline and
follow-up.

200 per arm. Baseline ~ normal(50, 10); follow-up = 20 + 0.6 x baseline - 3 x treated + noise (SD 6).
The true treatment effect on follow-up, adjusted for baseline, is -3, and the baseline coefficient
is 0.6.

The baseline coefficient is 0.6, neither 0 nor 1, on purpose. Analysing follow-up alone ignores the
baseline (right when the coefficient is 0); analysing the change score forces a coefficient of 1.
ANCOVA estimates it, which is why it is the most precise of the three in a trial. All three are
unbiased here because allocation is randomized; they differ in precision, which the controls show.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import random

SEED = 20261007
PER_ARM = 200
EFFECT = -3.0
BASELINE_COEF = 0.6


def rows(seed=SEED, per_arm=PER_ARM):
    rng = random.Random(seed)
    i = 0
    for treated in (1, 0):
        for _ in range(per_arm):
            i += 1
            baseline = rng.gauss(50.0, 10.0)
            followup = 20 + BASELINE_COEF * baseline + EFFECT * treated + rng.gauss(0.0, 6.0)
            yield {"id": i, "treated": treated, "baseline": round(baseline, 2), "followup": round(followup, 2)}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "treated", "baseline", "followup"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} rows")
