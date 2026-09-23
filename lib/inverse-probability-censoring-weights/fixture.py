"""Seeded fixture for inverse probability of censoring weights (IPCW): informative dropout in a trial.

2,000 patients randomized 1:1 (A). A baseline prognostic factor L ~ N(0, 1) raises the risk of the
outcome, logit P(Y = 1) = -1 + 1.5 L - 0.7 A, and also the chance of dropping out before the outcome is
measured, logit P(dropout) = -1.2 + 1.5 L - 1.5 A: sicker patients drop out more, and more so on
control, without the treatment's benefit. The outcome of the patients who stay is therefore not
representative of their arm, and the complete-case risk difference is biased (by about +0.05). Weighting each patient who stayed by 1 / P(stayed | A, L)
restores the arm (missing at random given A and L).

The target is the marginal risk difference E[Y(1)] - E[Y(0)], integrated from the data-generating
process by `truth()`. The fixture writes the outcome as blank for patients who dropped out.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261118
N = 2000


def expit(x):
    return 1.0 / (1.0 + math.exp(-x))


def risk(a, l):
    return expit(-1.0 + 1.5 * l - 0.7 * a)


def rows(seed=SEED):
    rng = random.Random(seed)
    arms = [1] * (N // 2) + [0] * (N // 2)
    rng.shuffle(arms)
    for i, a in enumerate(arms, start=1):
        l = rng.gauss(0.0, 1.0)
        y = 1 if rng.random() < risk(a, l) else 0
        dropped = rng.random() < expit(-1.2 + 1.5 * l - 1.5 * a)
        yield {"id": i, "treated": a, "l": round(l, 6), "event": "" if dropped else y}


def _integrate(f, lo=-10.0, hi=10.0, steps=20000):
    h = (hi - lo) / steps
    total = 0.0
    for k in range(steps + 1):
        x = lo + k * h
        w = 1 if k in (0, steps) else (4 if k % 2 else 2)
        total += w * f(x) * math.exp(-0.5 * x * x) / math.sqrt(2 * math.pi)
    return total * h / 3


def truth():
    r1, r0 = _integrate(lambda l: risk(1, l)), _integrate(lambda l: risk(0, l))
    return {"risk_treated": r1, "risk_control": r0, "risk_difference": r1 - r0}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "treated", "l", "event"])
        writer.writeheader()
        writer.writerows(data)
    t = truth()
    print(f"wrote fixture.csv: {len(data)} patients, {sum(1 for r in data if r['event'] == '')} dropped out; "
          f"risk difference {t['risk_difference']:.10f}")
