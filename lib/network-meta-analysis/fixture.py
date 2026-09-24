"""Seeded fixture for a network meta-analysis of three treatments from two-arm trials.

Placebo (A) and two active treatments (B, C). 15 two-arm trials: 6 compare A with B, 5 A with C, and 4
B with C. True log odds ratios against placebo: B -0.5, C -0.8, so C against B is -0.3 (the network is
consistent). Each trial's effect varies around these with between-trial SD 0.15. Arms have 100 to 300
patients and the placebo-arm risk is 0.3 (trial-specific baseline 0.25 to 0.35).

The data are one row per trial contrast: the log odds ratio of treat1 against treat2 and its Woolf
standard error.

The seed is 20261205, not the first one tried: at 20261204 the estimated tau^2 was 0 (Q below its df),
so the committed fixture would not have exercised the tau^2 estimator in the cross-language check. The
estimator is unbiased over seeds either way (control file); tau^2 is positive in about 77% of them.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261205
D = {"A": 0.0, "B": -0.5, "C": -0.8}
DESIGNS = [("B", "A")] * 6 + [("C", "A")] * 5 + [("C", "B")] * 4
TAU = 0.15


def expit(x):
    return 1.0 / (1.0 + math.exp(-x))


def truth():
    return {"log_or_B_vs_A": D["B"], "log_or_C_vs_A": D["C"], "log_or_C_vs_B": D["C"] - D["B"]}


def rows(seed=SEED):
    rng = random.Random(seed)
    for s, (t1, t2) in enumerate(DESIGNS, start=1):
        base = math.log(0.3 / 0.7) + rng.uniform(-0.25, 0.25)
        delta = D[t1] - D[t2] + rng.gauss(0.0, TAU)
        n = rng.randint(100, 300)
        p2 = expit(base + D[t2])
        p1 = expit(base + D[t2] + delta)
        e1 = sum(1 for _ in range(n) if rng.random() < p1)
        e2 = sum(1 for _ in range(n) if rng.random() < p2)
        te = math.log(e1 * (n - e2) / ((n - e1) * e2))
        se = math.sqrt(1 / e1 + 1 / (n - e1) + 1 / e2 + 1 / (n - e2))
        yield {"study": s, "treat1": t1, "treat2": t2, "TE": round(te, 10), "seTE": round(se, 10)}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["study", "treat1", "treat2", "TE", "seTE"])
        writer.writeheader()
        writer.writerows(data)
    print(f"wrote fixture.csv: {len(data)} two-arm trials; " + ", ".join(f"{k} {v}" for k, v in truth().items()))
