"""Seeded fixture for developing a penalized prediction model and correcting its apparent performance.

A development cohort of 300 patients with 10 candidate predictors X1..X10 ~ N(0, 1), of which four
matter: logit P(Y = 1) = -1.2 + 0.9 X1 + 0.6 X2 + 0.4 X3 + 0.2 X4. A model developed and assessed on the
same 300 patients looks better than it will in new ones (optimism). An external validation cohort of
5,000 patients from the same population, in the same file (cohort = "validation"), shows how the model
actually performs; the bootstrap-corrected estimates from the development data alone should predict it.

The development rows carry a fold number (1 to 5) for cross-validation, and the fixture supplies the 50
bootstrap resamples as row indices, so both languages use the same resamples; `bootstrap(seed)` returns
them and running the file writes bootstrap.csv beside fixture.csv.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSVs.
"""

import csv
import math
import random

SEED = 20261127
N_DEV = 300
N_VAL = 5000
P = 10
BETA = (-1.2, 0.9, 0.6, 0.4, 0.2)
FOLDS = 5
B = 50


def expit(x):
    return 1.0 / (1.0 + math.exp(-x))


def rows(seed=SEED):
    rng = random.Random(seed)
    for cohort, n in (("development", N_DEV), ("validation", N_VAL)):
        for i in range(n):
            x = [rng.gauss(0.0, 1.0) for _ in range(P)]
            eta = BETA[0] + sum(b * v for b, v in zip(BETA[1:], x))
            row = {"cohort": cohort, "fold": (i % FOLDS) + 1 if cohort == "development" else 0,
                   "event": 1 if rng.random() < expit(eta) else 0}
            row.update({f"x{j + 1}": round(v, 6) for j, v in enumerate(x)})
            yield row


def bootstrap(seed=SEED):
    """B resamples of the development rows, as 0-based indices, one list per resample."""
    rng = random.Random(seed + 1)
    return [[rng.randrange(N_DEV) for _ in range(N_DEV)] for _ in range(B)]


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(data[0]))
        writer.writeheader()
        writer.writerows(data)
    with open("bootstrap.csv", "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["resample", "rows"])
        for b, idx in enumerate(bootstrap(), start=1):
            writer.writerow([b, " ".join(map(str, idx))])
    dev = [r for r in data if r["cohort"] == "development"]
    print(f"wrote fixture.csv ({len(dev)} development, {len(data) - len(dev)} validation; "
          f"{sum(r['event'] for r in dev)} development events) and bootstrap.csv ({B} resamples)")
