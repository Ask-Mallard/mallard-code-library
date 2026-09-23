"""Seeded fixture for recurrent events: repeated hospital admissions over follow-up.

1000 people, treated with probability 0.5, followed for 1 to 3 years. Each person has their own
admission rate: 0.8 a year x a gamma frailty (mean 1, variance 0.5) x 0.7 if treated. Admissions then
arrive as a Poisson process at that rate, so one person can have several. Because the frailty is
independent of treatment, the population-average (marginal) rate ratio is 0.7, which is what the
Andersen-Gill model with a robust variance estimates.

Written in counting-process form: one row per gap between admissions, (start, stop], event = 1 if the
gap ended in an admission. Follow-up time after the last admission is a final row with event 0.

Some people are admitted far more often than others, so a person's rows are correlated. The model-
based (naive) standard error ignores that and is too small; the robust (sandwich) error clusters on
the person. Analysing only the FIRST admission discards most of the events.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261020
N = 1000
BASE_RATE = 0.8
RATE_RATIO = 0.7
LOG_RATE_RATIO = math.log(RATE_RATIO)
FRAILTY_SHAPE = 2.0   # gamma frailty with mean 1 and variance 1 / shape


def rows(seed=SEED, n=N):
    rng = random.Random(seed)
    for i in range(1, n + 1):
        treated = 1 if rng.random() < 0.5 else 0
        end = rng.uniform(1.0, 3.0)
        rate = BASE_RATE * rng.gammavariate(FRAILTY_SHAPE, 1 / FRAILTY_SHAPE) * (RATE_RATIO if treated else 1.0)
        start, t = 0.0, rng.expovariate(rate)
        while t < end:
            yield {"id": i, "treated": treated, "start": round(start, 7), "stop": round(t, 7), "event": 1}
            start, t = t, t + rng.expovariate(rate)
        yield {"id": i, "treated": treated, "start": round(start, 7), "stop": round(end, 7), "event": 0}


if __name__ == "__main__":
    data = [r for r in rows() if r["stop"] > r["start"]]
    assert len(data) == len(list(rows())), "a zero-length interval appeared; change the rounding"
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "treated", "start", "stop", "event"])
        writer.writeheader()
        writer.writerows(data)
    events = sum(r["event"] for r in data)
    print(f"wrote fixture.csv: {len(data)} rows, {events} admissions among {N} people")
