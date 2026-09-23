"""Seeded fixture for a self-controlled case series (SCCS): a vaccine and an acute adverse event.

Children are followed from day 0 to day 365 of the second year of life. Each is vaccinated on a day
between 120 and 240; the 28 days from vaccination are the RISK window. The event rate also varies with
age: days 0-119 x 1, days 120-239 x 2, days 240-364 x 0.5, so vaccination falls in the high-risk age
band and an analysis that ignores age credits the age effect to the vaccine. Each child has their own
underlying rate (gamma frailty, variance 0.5), which the design removes.

Events in each interval are Poisson(0.0015 x frailty x age multiplier x RR^risk x days), with RR = 3.
Only children with at least one event are kept, as an SCCS samples cases. The data are one row per
child per interval (age band x risk status), with its length in days and its event count.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261112
CASES = 400
BASE = 0.0015
AGE_BANDS = ((0, 120, 1.0), (120, 240, 2.0), (240, 365, 0.5))
RISK_DAYS = 28
LOG_RR = math.log(3.0)


def _poisson(rng, mean):
    # Knuth's method; the means here are small.
    limit, k, p = math.exp(-mean), 0, 1.0
    while True:
        p *= rng.random()
        if p <= limit:
            return k
        k += 1


def rows(seed=SEED):
    rng = random.Random(seed)
    child = 0
    while child < CASES:
        frailty = rng.gammavariate(2.0, 0.5)  # mean 1, variance 0.5
        vaccinated = rng.randint(120, 240)
        cuts = sorted({0, 365, vaccinated, min(vaccinated + RISK_DAYS, 365)} | {b[0] for b in AGE_BANDS})
        intervals = []
        for lo, hi in zip(cuts, cuts[1:]):
            band = next(i for i, b in enumerate(AGE_BANDS) if b[0] <= lo < b[1])
            risk = 1 if vaccinated <= lo < vaccinated + RISK_DAYS else 0
            mean = BASE * frailty * AGE_BANDS[band][2] * math.exp(LOG_RR * risk) * (hi - lo)
            intervals.append((band + 1, risk, hi - lo, _poisson(rng, mean)))
        if sum(iv[3] for iv in intervals) == 0:
            continue  # not a case: never sampled
        child += 1
        for band, risk, days, events in intervals:
            yield {"child": child, "age_band": band, "risk": risk, "days": days, "events": events}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["child", "age_band", "risk", "days", "events"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {CASES} cases, {len(data)} intervals, {sum(r['events'] for r in data)} events; "
          f"log RR {LOG_RR:.10f}")
