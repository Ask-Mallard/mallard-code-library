"""Seeded fixture for a Fine-Gray competing-risks regression, generated the way Fine and Gray (1999)
simulated their own model, so the subdistribution hazard ratio is known exactly.

1000 people, treated with probability 0.5. The cumulative incidence of cause 1 (relapse) is

    F1(t | x) = 1 - [1 - p (1 - exp(-t))] ^ exp(beta x),   p = 0.4, beta = -0.5,

which is exactly a proportional subdistribution-hazards model with log subdistribution hazard ratio
-0.5 for treatment. A person fails from cause 1 with probability F1(inf | x) = 1 - (1 - p)^exp(beta x);
their time then comes from F1 rescaled to that probability. Otherwise they fail from cause 2 (death) at
an exponential time with rate 1. Censoring is uniform on 0 to 4 years.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261019
N = 1000
P = 0.4
LOG_SHR = -0.5


def rows(seed=SEED, n=N):
    rng = random.Random(seed)
    for i in range(1, n + 1):
        treated = 1 if rng.random() < 0.5 else 0
        e = math.exp(LOG_SHR * treated)
        p1 = 1 - (1 - P) ** e
        if rng.random() < p1:
            u = rng.random()
            t = -math.log(1 - (1 - (1 - u * p1) ** (1 / e)) / P)
            cause = 1
        else:
            t = rng.expovariate(1.0)
            cause = 2
        censor = rng.uniform(0.0, 4.0)
        yield {"id": i, "treated": treated, "time": round(min(t, censor), 7), "status": cause if t <= censor else 0}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "treated", "time", "status"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} rows; cause 1 {sum(r['status'] == 1 for r in data)}, "
          f"cause 2 {sum(r['status'] == 2 for r in data)}, censored {sum(r['status'] == 0 for r in data)}")
