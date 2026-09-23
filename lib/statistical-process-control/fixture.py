"""Seeded fixture for statistical process control charts in a hospital quality programme.

24 months of data from one hospital:
- central-line infections over catheter-days (3,000 to 5,000 a month), at 2 per 1,000 catheter-days for
  months 1 to 18 and 5 per 1,000 from month 19 (a real change: a u-chart should signal);
- 30-day readmissions among discharges (150 to 250 a month), steady at 0.12 (a p-chart should not);
- mean length of stay, steady at 5.0 days with month-to-month SD 0.3 (an individuals chart should not).

Limits are set from the first 12 months (the baseline) and applied to all 24.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261130
MONTHS = 24
BASE_RATE, SHIFT_RATE, SHIFT_FROM = 2.0 / 1000, 5.0 / 1000, 19
P_READMIT = 0.12
MEAN_LOS, SD_LOS = 5.0, 0.3


def _poisson(rng, mean):
    # Inversion; the monthly means here are at most about 25, far from exp() underflow.
    k, p, cum, u = 0, math.exp(-mean), math.exp(-mean), rng.random()
    while u > cum:
        k += 1
        p *= mean / k
        cum += p
    return k


def rows(seed=SEED):
    rng = random.Random(seed)
    for m in range(1, MONTHS + 1):
        days = rng.randint(3000, 5000)
        rate = SHIFT_RATE if m >= SHIFT_FROM else BASE_RATE
        discharges = rng.randint(150, 250)
        readmissions = sum(1 for _ in range(discharges) if rng.random() < P_READMIT)
        yield {"month": m, "catheter_days": days, "infections": _poisson(rng, rate * days),
               "discharges": discharges, "readmissions": readmissions,
               "mean_los": round(rng.gauss(MEAN_LOS, SD_LOS), 6)}


def truth():
    return {"u_centre_per_1000": BASE_RATE * 1000, "p_centre": P_READMIT, "i_centre": MEAN_LOS}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(data[0]))
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} months; infections {sum(r['infections'] for r in data)}")
