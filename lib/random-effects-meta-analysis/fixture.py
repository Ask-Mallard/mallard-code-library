"""Seeded fixture for a random-effects meta-analysis of odds ratios from 15 trials.

Each trial randomizes between 60 and 300 patients per arm. The control-arm risk varies between trials
(0.2 to 0.4), and the true log odds ratio varies around -0.4 with between-trial SD 0.2 (tau^2 = 0.04):
trial i has log OR theta_i ~ N(-0.4, 0.2^2). Events in each arm are binomial. The summary estimand is
the MEAN of the distribution of true effects, -0.4.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV (one row per
trial, 2x2 counts).
"""

import csv
import math
import random

SEED = 20261202
TRIALS = 15
MU, TAU = -0.4, 0.2


def expit(x):
    return 1.0 / (1.0 + math.exp(-x))


def rows(seed=SEED):
    rng = random.Random(seed)
    for t in range(1, TRIALS + 1):
        n = rng.randint(60, 300)
        p0 = rng.uniform(0.2, 0.4)
        theta = rng.gauss(MU, TAU)
        p1 = expit(math.log(p0 / (1 - p0)) + theta)
        e1 = sum(1 for _ in range(n) if rng.random() < p1)
        e0 = sum(1 for _ in range(n) if rng.random() < p0)
        yield {"trial": t, "events_treated": e1, "n_treated": n, "events_control": e0, "n_control": n}


if __name__ == "__main__":
    data = list(rows())
    assert all(0 < r["events_treated"] < r["n_treated"] and 0 < r["events_control"] < r["n_control"] for r in data)
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(data[0]))
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {TRIALS} trials; mean log OR {MU}, tau {TAU}")
