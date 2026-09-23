"""Seeded fixture for an incidence rate: events over person-time in a cohort.

Each person is followed for a different length of time, uniform between 0.5 and 5 years, and
events arrive as a Poisson process at 0.08 per person-year. Events are counted by simulating the
exponential gaps between them, so the count is exactly Poisson given follow-up and the rate is a
property of the generator.

Follow-up varies on purpose. With equal follow-up a rate is a count divided by a constant, and an
implementation that divided events by the number of PEOPLE rather than by PERSON-TIME would still
land near the truth after rescaling. With unequal follow-up the two separate.

The rate is the same for everyone, so the Poisson assumption behind the exact interval holds here.
Real cohorts are often overdispersed (some people are prone to repeated events), and then the exact
Poisson interval is too narrow; README.md says what to use instead.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import random

SEED = 20260924
N = 600
RATE = 0.08  # events per person-year


def rows(seed=SEED, n=N, rate=RATE):
    rng = random.Random(seed)
    for i in range(1, n + 1):
        followup = round(rng.uniform(0.5, 5.0), 4)
        events, t = 0, rng.expovariate(rate)
        while t <= followup:
            events += 1
            t += rng.expovariate(rate)
        yield {"id": i, "followup_years": followup, "events": events}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "followup_years", "events"])
        writer.writeheader()
        writer.writerows(data)
    events = sum(r["events"] for r in data)
    pt = sum(r["followup_years"] for r in data)
    print(f"wrote fixture.csv: {len(data)} people, {events} events over {pt:.1f} person-years, "
          f"observed rate {events / pt:.4f} (true {RATE})")
