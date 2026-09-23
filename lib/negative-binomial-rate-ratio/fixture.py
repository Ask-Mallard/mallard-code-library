"""Seeded fixture for overdispersed event counts over unequal follow-up.

1000 people, treated with probability 0.5, followed for 0.5 to 3 years. Each person has their own
event rate, drawn from a gamma distribution around 0.8 per year (control) or 0.8 x 0.7 (treated) with
shape theta = 2; events then arrive as a Poisson process at that person's rate. That gamma-Poisson
mixture IS the negative binomial (NB2) model, with true rate ratio 0.7 (log -0.357) and theta 2.

Some people are prone to events and some are not, which is the usual reason clinical counts (falls,
admissions, exacerbations) are overdispersed.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261009
N = 1000
RATE_CONTROL = 0.8
RATE_RATIO = 0.7
THETA = 2.0
LOG_RATE_RATIO = math.log(RATE_RATIO)


def rows(seed=SEED, n=N):
    rng = random.Random(seed)
    for i in range(1, n + 1):
        treated = 1 if rng.random() < 0.5 else 0
        followup = round(rng.uniform(0.5, 3.0), 4)
        mean = RATE_CONTROL * (RATE_RATIO if treated else 1.0) * followup
        personal = rng.gammavariate(THETA, mean / THETA)   # gamma mean = mean
        events, t = 0, (rng.expovariate(personal) if personal > 0 else math.inf)
        while t <= 1.0:
            events += 1
            t += rng.expovariate(personal)
        yield {"id": i, "treated": treated, "followup_years": followup, "events": events}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "treated", "followup_years", "events"])
        writer.writeheader()
        writer.writerows(data)
    ev = [r["events"] for r in data]
    m = sum(ev) / len(ev)
    v = sum((e - m) ** 2 for e in ev) / (len(ev) - 1)
    print(f"wrote fixture.csv: {len(data)} rows, {sum(ev)} events, count mean {m:.2f} variance {v:.2f}")
