"""Seeded fixture for restricted mean survival time (RMST) in a two-arm trial.

700 per arm; survival is exponential with hazard 0.20 a year (control) and 0.12 (treated), and
follow-up ends administratively between 3 and 6 years. For an exponential, the mean survival up to a
horizon tau is (1 - exp(-hazard x tau)) / hazard, so at tau = 4 years the true RMST is 2.7534 years
(control) and 3.1768 (treated): the treated live 0.4234 years longer on average over those 4 years.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261017
PER_ARM = 700
HAZARD = {0: 0.20, 1: 0.12}
TAU = 4.0


def rmst(hazard, tau=TAU):
    return (1 - math.exp(-hazard * tau)) / hazard


TRUE_RMST = {arm: rmst(h) for arm, h in HAZARD.items()}
TRUE_DIFFERENCE = TRUE_RMST[1] - TRUE_RMST[0]


def rows(seed=SEED, per_arm=PER_ARM):
    rng = random.Random(seed)
    i = 0
    for arm in (1, 0):
        for _ in range(per_arm):
            i += 1
            t = rng.expovariate(HAZARD[arm])
            censor = rng.uniform(3.0, 6.0)
            yield {"id": i, "arm": arm, "time": round(min(t, censor), 5), "event": 1 if t <= censor else 0}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "arm", "time", "event"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} rows, {sum(r['event'] for r in data)} events; "
          f"true RMST {TRUE_RMST[1]:.4f} vs {TRUE_RMST[0]:.4f}, difference {TRUE_DIFFERENCE:.4f}")
