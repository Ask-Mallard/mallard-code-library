"""Seeded fixture for a case-crossover study: a transient exposure and an acute event.

Each of 600 patients had an acute event (say, a myocardial infarction). For each, exposure (say, heavy
exertion) is recorded in the hazard window just before the event and in 3 referent windows at other
times in the same month (time-stratified design). Each person's usual exposure probability differs,
p ~ Beta(2, 5), and is the kind of between-person confounder the design removes by comparing each
person with themselves.

The four windows' exposures are Bernoulli(p); the event falls in one of the four windows with
probability proportional to exp(beta x), beta = log(2.5), which is the case-crossover likelihood. The
conditional odds ratio for exposure is 2.5.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261111
N = 600
WINDOWS = 4
LOG_OR = math.log(2.5)


def rows(seed=SEED):
    rng = random.Random(seed)
    for i in range(1, N + 1):
        p = rng.betavariate(2, 5)
        x = [1 if rng.random() < p else 0 for _ in range(WINDOWS)]
        weights = [math.exp(LOG_OR * v) for v in x]
        u = rng.random() * sum(weights)
        hazard = 0
        while u > weights[hazard]:
            u -= weights[hazard]
            hazard += 1
        for w in range(WINDOWS):
            yield {"patient": i, "window": w + 1, "hazard_window": 1 if w == hazard else 0, "exposed": x[w]}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["patient", "window", "hazard_window", "exposed"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {N} patients x {WINDOWS} windows; log OR {LOG_OR:.10f}")
