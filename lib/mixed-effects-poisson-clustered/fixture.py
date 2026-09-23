"""Seeded fixture for clustered counts: infections per patient-days in 40 wards.

40 wards, 20 per arm, 40 patients each with a length of stay of 2 to 20 days. Each ward has its own
infection rate: 0.03 a patient-day x a ward multiplier exp(N(0, 0.25)) x 0.6 if the ward is in the
intervention arm. Infections are Poisson given the ward's rate and the patient's days.

With a log link and a normal random intercept, the ward-conditional rate ratio (0.6) is also the
population-average one: the Poisson rate ratio is collapsible over a random intercept, unlike the odds
ratio of mixed-effects-logistic-clustered.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261027
WARDS_PER_ARM = 20
PER_WARD = 40
BASE_RATE = 0.03
RATE_RATIO = 0.6
LOG_RATE_RATIO = math.log(RATE_RATIO)
WARD_SD = 0.25


def rows(seed=SEED):
    rng = random.Random(seed)
    ward = 0
    for arm in (1, 0):
        for _ in range(WARDS_PER_ARM):
            ward += 1
            rate = BASE_RATE * math.exp(rng.gauss(0.0, WARD_SD)) * (RATE_RATIO if arm else 1.0)
            for _ in range(PER_WARD):
                days = rng.randint(2, 20)
                count, t = 0, rng.expovariate(rate)
                while t <= days:
                    count += 1
                    t += rng.expovariate(rate)
                yield {"ward": ward, "intervention": arm, "days": days, "infections": count}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["ward", "intervention", "days", "infections"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} patients in {2 * WARDS_PER_ARM} wards, "
          f"{sum(r['infections'] for r in data)} infections over {sum(r['days'] for r in data)} patient-days")
