"""Seeded fixture for multiplicity: a trial analysed on ten secondary outcomes.

200 patients randomized 1:1; ten continuous secondary outcomes, each with SD 1. The treatment shifts
the first three by 0.5, 0.4 and 0.3 SD; the other seven are unaffected. Each outcome is compared by
Welch's t-test, and the ten p-values are then adjusted. With seven true nulls, testing each at 0.05
unadjusted rejects at least one null with probability about 1 - 0.95^7 = 0.30.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import random

SEED = 20261129
N = 200
EFFECTS = (0.5, 0.4, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)


def rows(seed=SEED):
    rng = random.Random(seed)
    arms = [1] * (N // 2) + [0] * (N // 2)
    rng.shuffle(arms)
    for i, a in enumerate(arms, start=1):
        row = {"id": i, "treated": a}
        for j, e in enumerate(EFFECTS, start=1):
            row[f"y{j}"] = round(e * a + rng.gauss(0.0, 1.0), 6)
        yield row


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(data[0]))
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} patients, {len(EFFECTS)} outcomes, {sum(e > 0 for e in EFFECTS)} with an effect")
