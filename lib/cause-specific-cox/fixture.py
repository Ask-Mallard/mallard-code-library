"""Seeded fixture for competing risks: relapse (cause 1) and death without relapse (cause 2).

2500 people, treated with probability 0.5, followed up to 5 years with random censoring. Cause-specific
hazards are constant: relapse 0.10 a year (x 0.7 if treated) and death 0.06 a year (no treatment
effect). So the true cause-specific hazard ratio for relapse is 0.7.

A cause-specific hazard ratio is about the RATE of relapse among those still event-free. It is not the
effect on the cumulative incidence (the proportion who relapse by a time point), which also depends on
the competing death rate: see competing-risks-cumulative-incidence and the Fine-Gray entry for that
question.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261016
N = 2500
RELAPSE_HAZARD = 0.10
DEATH_HAZARD = 0.06
CAUSE_HR = 0.7
LOG_CS_HR = math.log(CAUSE_HR)


def rows(seed=SEED, n=N):
    rng = random.Random(seed)
    for i in range(1, n + 1):
        treated = 1 if rng.random() < 0.5 else 0
        h1 = RELAPSE_HAZARD * (CAUSE_HR if treated else 1.0)
        t = rng.expovariate(h1 + DEATH_HAZARD)
        cause = 1 if rng.random() < h1 / (h1 + DEATH_HAZARD) else 2
        censor = min(5.0, rng.expovariate(0.05))
        status = cause if t <= censor else 0
        yield {"id": i, "treated": treated, "time": round(min(t, censor), 7), "status": status}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "treated", "time", "status"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} rows; relapses {sum(r['status'] == 1 for r in data)}, "
          f"deaths {sum(r['status'] == 2 for r in data)}, censored {sum(r['status'] == 0 for r in data)}")
