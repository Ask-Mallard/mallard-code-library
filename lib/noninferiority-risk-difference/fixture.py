"""Seeded fixture for a noninferiority trial with a binary outcome and protocol deviations.

800 patients randomized 1:1 to a new treatment or the standard. Success with the standard is 0.80 and
with the new treatment 0.78, so the true risk difference (new minus standard) is -0.02. The prespecified
noninferiority margin is -0.10: the new treatment is noninferior if the lower 95% limit of the difference
lies above -0.10.

10% of the patients randomized to the new treatment receive the standard instead (their success is
0.80). An intention-to-treat analysis counts them as randomized, which pulls the difference towards 0
(-0.018) and so towards a noninferiority conclusion; the per-protocol analysis excludes them (-0.02).
Both are reported, as CONSORT's noninferiority extension asks.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import random

SEED = 20261128
N = 800
P_STANDARD, P_NEW = 0.80, 0.78
SWITCH = 0.10
MARGIN = -0.10


def truth():
    itt = (1 - SWITCH) * P_NEW + SWITCH * P_STANDARD - P_STANDARD
    return {"risk_difference_itt": itt, "risk_difference_pp": P_NEW - P_STANDARD}


def rows(seed=SEED):
    rng = random.Random(seed)
    arms = [1] * (N // 2) + [0] * (N // 2)
    rng.shuffle(arms)
    for i, a in enumerate(arms, start=1):
        switched = 1 if a == 1 and rng.random() < SWITCH else 0
        p = P_NEW if a == 1 and not switched else P_STANDARD
        yield {"id": i, "randomized_new": a, "per_protocol": 0 if switched else 1, "success": 1 if rng.random() < p else 0}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "randomized_new", "per_protocol", "success"])
        writer.writeheader()
        writer.writerows(data)
    t = truth()
    print(f"wrote fixture.csv: {len(data)} patients; ITT {t['risk_difference_itt']:.4f}, PP {t['risk_difference_pp']:.4f}, margin {MARGIN}")
