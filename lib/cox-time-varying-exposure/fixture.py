"""Seeded fixture for an exposure that starts during follow-up (for example, a drug started later).

1500 people followed for up to 5 years (plus random censoring). Each may start the exposure at a random
time (exponential, rate 0.3 a year; many never do within follow-up). The event hazard is 0.15 a year
while unexposed and 0.15 x 0.6 after starting, so the true hazard ratio for CURRENT exposure is 0.6.

Written in counting-process form: one row per person-period, (start, stop], with the exposure as it was
during that period. A person who starts the exposure contributes an unexposed row and then an exposed
one.

Coding exposure as "ever exposed" at baseline credits the time BEFORE the start to the exposed group:
time during which, by construction, they could not have had an exposed event. That is immortal time
bias, and it makes the exposure look protective whatever its true effect.
harness/survival_examples.test.py shows it.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261015
N = 1500
BASE_HAZARD = 0.15
HAZARD_RATIO = 0.6
LOG_HR = math.log(HAZARD_RATIO)
START_RATE = 0.3
ADMIN_END = 5.0


def rows(seed=SEED, n=N):
    rng = random.Random(seed)
    for i in range(1, n + 1):
        start = rng.expovariate(START_RATE)
        censor = min(ADMIN_END, rng.expovariate(0.05))
        unexposed_event = rng.expovariate(BASE_HAZARD)
        if unexposed_event < start:
            event_time = unexposed_event
        else:
            event_time = start + rng.expovariate(BASE_HAZARD * HAZARD_RATIO)
        end = round(min(event_time, censor), 7)
        event = 1 if event_time <= censor else 0
        s = round(start, 7)
        if 0 < s < end:
            yield {"id": i, "start": 0.0, "stop": s, "exposed": 0, "event": 0}
            yield {"id": i, "start": s, "stop": end, "exposed": 1, "event": event}
        else:
            yield {"id": i, "start": 0.0, "stop": end, "exposed": 0, "event": event}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "start", "stop", "exposed", "event"])
        writer.writeheader()
        writer.writerows(data)
    people = len({r["id"] for r in data})
    started = sum(1 for r in data if r["exposed"] == 1)
    print(f"wrote fixture.csv: {len(data)} rows for {people} people, {started} started the exposure, "
          f"{sum(r['event'] for r in data)} events")
