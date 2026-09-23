"""Seeded fixture for a method-comparison study: Bland-Altman limits of agreement.

150 patients measured by a reference method A and a new point-of-care method B. True values
T ~ N(100, 15^2); A = T + N(0, 4^2); B = T + 2 + N(0, 3^2). So B reads 2 units higher on average (bias
2), the differences B - A have SD 5, and 95% of differences lie within 2 +/- 1.96 x 5, i.e. -7.8 to 11.8.
There is no proportional bias. The two methods correlate at about 0.95, which says nothing about whether
they agree.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261125
N = 150
BIAS = 2.0
SD_A, SD_B = 4.0, 3.0


def truth():
    sd = math.hypot(SD_A, SD_B)
    return {"bias": BIAS, "sd_difference": sd, "loa_lower": BIAS - 1.959963984540054 * sd,
            "loa_upper": BIAS + 1.959963984540054 * sd}


def rows(seed=SEED):
    rng = random.Random(seed)
    for i in range(1, N + 1):
        t = rng.gauss(100.0, 15.0)
        yield {"id": i, "method_a": round(t + rng.gauss(0.0, SD_A), 6), "method_b": round(t + BIAS + rng.gauss(0.0, SD_B), 6)}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "method_a", "method_b"])
        writer.writeheader()
        writer.writerows(data)
    t = truth()
    print(f"wrote fixture.csv: {len(data)} patients; bias {t['bias']}, limits {t['loa_lower']:.6f} to {t['loa_upper']:.6f}")
