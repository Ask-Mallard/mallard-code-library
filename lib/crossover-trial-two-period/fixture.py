"""Seeded fixture for a two-period, two-treatment (AB/BA) crossover trial.

Each patient receives both treatments, in sequence AB or BA, with a washout between. 35 patients follow
AB and 25 follow BA: unequal on purpose. The outcome in each period is
patient level (SD 6) + 4 x period 2 + (-3) x treatment B + noise (SD 3), so the true treatment effect
(B minus A) is -3 and scores are 4 higher in period 2 whatever the treatment. No carryover.

A paired comparison of B with A that ignores the period is biased when sequences are unequal: the
period effect leaks into it in proportion to the imbalance (here 4 x (35 - 25) / 60 = 0.67).
The period-adjusted estimate compares the within-patient period differences between sequences.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import random

SEED = 20261024
SEQUENCES = {"AB": 35, "BA": 25}
PERIOD_EFFECT = 4.0
TREATMENT_EFFECT = -3.0   # B minus A


def rows(seed=SEED):
    rng = random.Random(seed)
    i = 0
    for seq, count in SEQUENCES.items():
        for _ in range(count):
            i += 1
            level = rng.gauss(50.0, 6.0)
            for period in (1, 2):
                treatment = seq[period - 1]
                y = level + PERIOD_EFFECT * (period == 2) + TREATMENT_EFFECT * (treatment == "B") + rng.gauss(0.0, 3.0)
                yield {"id": i, "sequence": seq, "period": period, "treatment": treatment, "y": round(y, 4)}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "sequence", "period", "treatment", "y"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} rows for {sum(SEQUENCES.values())} patients")
