"""Seeded fixture for bootstrap confidence intervals: the ratio of mean length of stay.

150 patients at each of two hospitals. Length of stay is right-skewed: log-normal with log-scale mean
1.6 at hospital A and 1.1 at hospital B, SD 0.6 at both. The ratio of MEANS (the quantity that drives
bed-days and costs) is exp(1.6 + 0.18) / exp(1.1 + 0.18) = exp(0.5) = 1.6487. With skewed data and a
ratio, the bootstrap is the usual route to an interval.

A median would be the wrong statistic for the BCa interval's jackknife acceleration: leaving one of 80
values out gives only two distinct medians, symmetrically, so the acceleration is exactly 0 (found when
this entry first used the difference in medians).

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261201
N_PER_GROUP = 150
MU_A, MU_B, SIGMA = 1.6, 1.1, 0.6


def truth():
    return {"mean_ratio": math.exp(MU_A - MU_B)}


def rows(seed=SEED):
    rng = random.Random(seed)
    for hospital, mu in (("A", MU_A), ("B", MU_B)):
        for i in range(N_PER_GROUP):
            yield {"hospital": hospital, "los": round(math.exp(rng.gauss(mu, SIGMA)), 6)}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["hospital", "los"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} patients; ratio of means {truth()['mean_ratio']:.10f}")
