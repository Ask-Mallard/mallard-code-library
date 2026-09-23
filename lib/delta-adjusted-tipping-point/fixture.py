"""Seeded fixture for a delta-adjusted tipping-point analysis: a trial with missing outcomes.

1,200 patients randomized 1:1. A baseline covariate X ~ N(0, 1) predicts the outcome:
Y = 10 + 1.5 A + 3 X + noise (SD 1.5). The outcome is missing for about a third of patients, more often
when X is high and more often in the treated arm: logit P(Y missing) = -1.5 + 1.5 X + 1.0 A. That is
missing at random given X and the arm, so an analysis that models Y on X within each arm is unbiased
(the effect is 1.5), while the complete-case difference in means is not.

Whether the data are MAR cannot be checked. The tipping-point analysis asks how much worse the treated
patients with missing outcomes would have to be, by a shift delta added to their predicted outcomes, for
the effect to lose statistical significance. The fixture writes y as blank where it is missing.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261116
N = 1200
EFFECT = 1.5


def expit(x):
    return 1.0 / (1.0 + math.exp(-x))


def rows(seed=SEED):
    rng = random.Random(seed)
    arms = [1] * (N // 2) + [0] * (N // 2)
    rng.shuffle(arms)
    for i, a in enumerate(arms, start=1):
        x = rng.gauss(0.0, 1.0)
        y = 10 + EFFECT * a + 3 * x + rng.gauss(0.0, 1.5)
        missing = rng.random() < expit(-1.5 + 1.5 * x + 1.0 * a)
        yield {"id": i, "treated": a, "x": round(x, 6), "y": "" if missing else round(y, 6)}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "treated", "x", "y"])
        writer.writeheader()
        writer.writerows(data)
    miss = sum(1 for r in data if r["y"] == "")
    print(f"wrote fixture.csv: {len(data)} patients, outcome missing for {miss}; effect {EFFECT}")
