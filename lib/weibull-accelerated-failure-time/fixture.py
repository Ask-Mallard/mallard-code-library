"""Seeded fixture for a parametric (Weibull) survival model.

600 people, treated with probability 0.5. Survival times are Weibull with shape 1.5 and scale 5 years
for controls; treatment multiplies every survival time by exp(0.4) (an accelerated failure time model),
so the true log time ratio is 0.4. For a Weibull, the same model is a proportional-hazards model with
log hazard ratio -0.4 x shape = -0.6. Follow-up ends uniformly between 2 and 8 years.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261018
N = 600
SHAPE = 1.5
SCALE = 5.0
LOG_TIME_RATIO = 0.4
LOG_HAZARD_RATIO = -LOG_TIME_RATIO * SHAPE


def rows(seed=SEED, n=N):
    rng = random.Random(seed)
    for i in range(1, n + 1):
        treated = 1 if rng.random() < 0.5 else 0
        u = rng.random()
        t = SCALE * math.exp(LOG_TIME_RATIO * treated) * (-math.log(u)) ** (1 / SHAPE)
        censor = rng.uniform(2.0, 8.0)
        yield {"id": i, "treated": treated, "time": round(min(t, censor), 5), "event": 1 if t <= censor else 0}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "treated", "time", "event"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} rows, {sum(r['event'] for r in data)} events")
