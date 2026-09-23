"""Seeded fixture for a case-cohort study: a random subcohort plus every case from a larger cohort.

A cohort of 4000 is followed up to 5 years; 30% are exposed and the event hazard is 0.03 a year x 1.8
if exposed, so the true log hazard ratio is log(1.8). Measuring the exposure in everyone is expensive,
so it is measured only in a random 10% SUBCOHORT and in every CASE. The committed CSV is that
case-cohort sample: subcohort members (cases or not) and cases outside the subcohort.

Cases are over-represented in the sample, so an ordinary Cox model on it treats the sample as if it
were the cohort. A case-cohort analysis weights the subcohort to stand for the whole cohort and uses a
robust variance.

`cohort(seed)` generates the FULL cohort, which a real study never sees; the controls use it to check
the case-cohort estimate against the full-cohort estimate. `rows(seed)` is the case-cohort sample.
Stdlib only. Running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261021
COHORT = 4000
SUBCOHORT_FRACTION = 0.10
BASE_HAZARD = 0.03
HAZARD_RATIO = 1.8
LOG_HR = math.log(HAZARD_RATIO)


def cohort(seed=SEED, n=COHORT):
    rng = random.Random(seed)
    subcohort = set(rng.sample(range(1, n + 1), int(SUBCOHORT_FRACTION * n)))
    for i in range(1, n + 1):
        exposed = 1 if rng.random() < 0.3 else 0
        t = rng.expovariate(BASE_HAZARD * (HAZARD_RATIO if exposed else 1.0))
        censor = min(5.0, rng.expovariate(0.02))
        yield {"id": i, "exposed": exposed, "time": round(min(t, censor), 7),
               "event": 1 if t <= censor else 0, "subcohort": 1 if i in subcohort else 0}


def rows(seed=SEED, n=COHORT):
    for r in cohort(seed, n):
        if r["subcohort"] or r["event"]:
            yield r


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "exposed", "time", "event", "subcohort"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} of {COHORT}; subcohort {sum(r['subcohort'] for r in data)}, "
          f"cases {sum(r['event'] for r in data)}")
