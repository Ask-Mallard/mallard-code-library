"""Seeded fixture for a controlled interrupted time series: an intervention site and a control site.

60 months at two sites; the intervention starts at month 37 at one of them. At the SAME month something
else changes at both sites (a coding change, a regional campaign): a level drop of 2 everywhere. Both
series also share month-to-month shocks.

intervention: y = 50 + 0.2 t - 2 [t >= 37] - 4 [t >= 37] - 0.3 (t - 36)+ + s_t + e_t
control:      y = 45 + 0.25 t - 2 [t >= 37]                             + s_t + f_t

s_t is a shared AR(1) shock (rho 0.5, SD 1); e_t and f_t are site-specific AR(1) noise (rho 0.3, SD 0.5).
The intervention's effect is a level change of -4 and a slope change of -0.3 a month RELATIVE to the
control. A single-series analysis of the intervention site credits the shared drop to the intervention
(-6); the control series removes it.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261108
MONTHS = 60
START = 37
LEVEL_CHANGE = -4.0
SLOPE_CHANGE = -0.3
SHARED_DROP = -2.0


def rows(seed=SEED):
    rng = random.Random(seed)
    s = rng.gauss(0.0, 1.0 / math.sqrt(1 - 0.25))
    e = rng.gauss(0.0, 0.5 / math.sqrt(1 - 0.09))
    f = rng.gauss(0.0, 0.5 / math.sqrt(1 - 0.09))
    for t in range(1, MONTHS + 1):
        if t > 1:
            s = 0.5 * s + rng.gauss(0.0, 1.0)
            e = 0.3 * e + rng.gauss(0.0, 0.5)
            f = 0.3 * f + rng.gauss(0.0, 0.5)
        post = 1 if t >= START else 0
        since = max(t - (START - 1), 0)
        yi = 50 + 0.2 * t + SHARED_DROP * post + LEVEL_CHANGE * post + SLOPE_CHANGE * since + s + e
        yc = 45 + 0.25 * t + SHARED_DROP * post + s + f
        yield {"month": t, "post": post, "months_since": since,
               "intervention_rate": round(yi, 6), "control_rate": round(yc, 6)}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["month", "post", "months_since", "intervention_rate", "control_rate"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {MONTHS} months x 2 sites; level {LEVEL_CHANGE}, slope {SLOPE_CHANGE}, shared drop {SHARED_DROP}")
