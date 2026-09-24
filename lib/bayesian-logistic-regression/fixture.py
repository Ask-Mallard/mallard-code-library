"""Seeded fixture for a Bayesian logistic regression with stated priors.

1,000 patients in a randomized trial; age (standardized) is a prognostic covariate.
logit P(Y = 1) = -1.0 + 0.6 x treated + 0.3 x age. The quantity of interest is the treatment's
conditional log odds ratio, 0.6, and the posterior probability that it is above 0.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261206
N = 1000
BETA = (-1.0, 0.6, 0.3)


def expit(x):
    return 1.0 / (1.0 + math.exp(-x))


def truth():
    return {"treatment_log_or": BETA[1]}


def rows(seed=SEED):
    rng = random.Random(seed)
    arms = [1] * (N // 2) + [0] * (N // 2)
    rng.shuffle(arms)
    for i, a in enumerate(arms, start=1):
        age = rng.gauss(0.0, 1.0)
        y = 1 if rng.random() < expit(BETA[0] + BETA[1] * a + BETA[2] * age) else 0
        yield {"id": i, "treated": a, "age": round(age, 6), "event": y}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "treated", "age", "event"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} patients, {sum(r['event'] for r in data)} events; treatment log OR {BETA[1]}")
