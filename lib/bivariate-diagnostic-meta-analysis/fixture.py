"""Seeded fixture for a meta-analysis of diagnostic accuracy studies: the bivariate (Reitsma) model.

20 studies of one test. Study i's true sensitivity and specificity come from a bivariate normal on the
logit scale: mean logit sensitivity 1.5 (0.818) and mean logit specificity 2.0 (0.881), between-study
SDs 0.5 and 0.6, correlation -0.4 between logit sensitivity and logit specificity (studies with a lower
threshold are more sensitive and less specific). Each study has 30 to 150 patients with the condition
and 100 to 400 without; TP and TN are binomial.

The model's parameters are the MEANS of the logit sensitivity and specificity across studies; `truth()`
returns them.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261203
STUDIES = 20
MU_SENS, MU_SPEC = 1.5, 2.0
SD_SENS, SD_SPEC, RHO = 0.5, 0.6, -0.4


def expit(x):
    return 1.0 / (1.0 + math.exp(-x))


def truth():
    return {"logit_sensitivity": MU_SENS, "logit_specificity": MU_SPEC}


def rows(seed=SEED):
    rng = random.Random(seed)
    for s in range(1, STUDIES + 1):
        z1, z2 = rng.gauss(0.0, 1.0), rng.gauss(0.0, 1.0)
        ls = MU_SENS + SD_SENS * z1
        lp = MU_SPEC + SD_SPEC * (RHO * z1 + math.sqrt(1 - RHO ** 2) * z2)
        nd, nn = rng.randint(30, 150), rng.randint(100, 400)
        tp = sum(1 for _ in range(nd) if rng.random() < expit(ls))
        tn = sum(1 for _ in range(nn) if rng.random() < expit(lp))
        yield {"study": s, "TP": tp, "FN": nd - tp, "FP": nn - tn, "TN": tn}


if __name__ == "__main__":
    data = list(rows())
    assert all(min(r["TP"], r["FN"], r["FP"], r["TN"]) > 0 for r in data), "a zero cell would invoke a continuity correction"
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["study", "TP", "FN", "FP", "TN"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {STUDIES} studies; logit sensitivity {MU_SENS}, logit specificity {MU_SPEC}")
