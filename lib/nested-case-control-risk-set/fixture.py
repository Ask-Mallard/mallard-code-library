"""Seeded fixture for a nested case-control study: controls sampled from each case's risk set.

A cohort of 3,000 is followed for up to 5 years. Exposure X ~ Bernoulli(0.3) and age (standardized,
N(0, 1)) set the hazard: lambda = 0.08 x exp(log(2) X + 0.5 age) a year, exponential event times,
administrative censoring at 5 years. The outcome is COMMON (about 40% have it within 5 years), on
purpose. The cohort itself is never analysed: for each case, 4 controls are drawn at random from those
still at risk at the case's event time (incidence density sampling), which is what a nested case-control
study does when measuring exposure on the whole cohort is too expensive.

Analysed with conditional logistic regression on the matched sets, the odds ratio estimates the HAZARD
RATIO of the full cohort, here 2, with no rare-disease assumption. A person can be a control for several
cases and later a case. `cohort(rng)` is exposed so the control file can draw controls the wrong way
(from everyone event-free at the end) and show that that odds ratio does not estimate the hazard ratio
when the outcome is common.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import bisect
import csv
import math
import random

SEED = 20261113
COHORT = 3000
CONTROLS = 4
BASE = 0.08
LOG_HR = math.log(2.0)
AGE_BETA = 0.5
FOLLOW_UP = 5.0


def cohort(rng):
    """(time, event, exposed, age) for every cohort member, drawn from rng."""
    people = []
    for _ in range(COHORT):
        x = 1 if rng.random() < 0.3 else 0
        age = rng.gauss(0.0, 1.0)
        t = rng.expovariate(BASE * math.exp(LOG_HR * x + AGE_BETA * age))
        people.append((min(t, FOLLOW_UP), t < FOLLOW_UP, x, round(age, 6)))
    return people


def rows(seed=SEED):
    rng = random.Random(seed)
    people = cohort(rng)
    order = sorted(range(COHORT), key=lambda i: people[i][0])
    times = [people[i][0] for i in order]
    case_set = 0
    for i in order:
        time, event, x, age = people[i]
        if not event:
            continue
        case_set += 1
        yield {"set": case_set, "case": 1, "person": i + 1, "exposed": x, "age": age}
        at_risk = order[bisect.bisect_right(times, time):]  # still under follow-up after this time
        for j in rng.sample(at_risk, CONTROLS):
            yield {"set": case_set, "case": 0, "person": j + 1, "exposed": people[j][2], "age": people[j][3]}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["set", "case", "person", "exposed", "age"])
        writer.writeheader()
        writer.writerows(data)
    sets = max(r["set"] for r in data)
    print(f"wrote fixture.csv: {sets} matched sets of 1 case and {CONTROLS} controls from a cohort of {COHORT}; "
          f"log HR {LOG_HR:.10f}")
