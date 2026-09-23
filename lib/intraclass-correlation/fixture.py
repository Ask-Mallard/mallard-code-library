"""Seeded fixture for inter-rater reliability of a continuous measurement: the intraclass correlation.

60 patients each measured once by the same 4 raters (for example, a lesion diameter in mm). The
measurement is y = 50 + subject effect (SD 8) + rater effect (SD 3) + error (SD 4). Raters differ
systematically (one reads larger than another), so ABSOLUTE AGREEMENT and CONSISTENCY differ:

ICC(A,1) = 64 / (64 + 9 + 16) = 0.7191   (two-way random, absolute agreement, single rater)
ICC(C,1) = 64 / (64 + 16)     = 0.8000   (consistency: rater offsets ignored)

`truth()` returns both. The rater effects are drawn per fixture, so an estimate of ICC(A,1) from 4
raters is noisier than one of ICC(C,1).

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV (wide: one column
per rater).
"""

import csv
import random

SEED = 20261124
SUBJECTS = 60
RATERS = 4
SD_SUBJECT, SD_RATER, SD_ERROR = 8.0, 3.0, 4.0


def truth():
    s, r, e = SD_SUBJECT ** 2, SD_RATER ** 2, SD_ERROR ** 2
    return {"icc_agreement": s / (s + r + e), "icc_consistency": s / (s + e)}


def rows(seed=SEED):
    rng = random.Random(seed)
    raters = [rng.gauss(0.0, SD_RATER) for _ in range(RATERS)]
    for i in range(1, SUBJECTS + 1):
        s = rng.gauss(0.0, SD_SUBJECT)
        row = {"subject": i}
        for j, r in enumerate(raters, start=1):
            row[f"rater_{j}"] = round(50 + s + r + rng.gauss(0.0, SD_ERROR), 6)
        yield row


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(data[0]))
        writer.writeheader()
        writer.writerows(data)
    t = truth()
    print(f"wrote fixture.csv: {SUBJECTS} subjects x {RATERS} raters; ICC(A,1) {t['icc_agreement']:.10f}, "
          f"ICC(C,1) {t['icc_consistency']:.10f}")
