"""Seeded fixture for sensitivity to unmeasured confounding: E-value and a simple bias analysis.

A cohort of 4,000. An UNMEASURED binary factor U ~ Bernoulli(0.4) makes exposure more likely,
logit P(A = 1) = -0.5 + 1.0 U, and raises the risk of the outcome. The risk follows a log-linear model
with no interaction: P(Y = 1) = 0.10 x 1.5^A x 2.5^U. The causal risk ratio for A is 1.5; the crude
risk ratio, which cannot adjust for U, is larger.

With this model the classic bias factor for a binary unmeasured confounder,
(p1 (RR_UY - 1) + 1) / (p0 (RR_UY - 1) + 1), with p1 = P(U = 1 | A = 1), p0 = P(U = 1 | A = 0) and
RR_UY = 2.5, is exact, so dividing the crude RR by it recovers 1.5 in expectation. `bias_parameters()`
returns those true values; an analyst would have to supply them from external evidence.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261117
N = 4000
P_U = 0.4
RR_A = 1.5
RR_UY = 2.5
BASE = 0.10


def expit(x):
    return 1.0 / (1.0 + math.exp(-x))


def rows(seed=SEED):
    rng = random.Random(seed)
    for i in range(1, N + 1):
        u = 1 if rng.random() < P_U else 0
        a = 1 if rng.random() < expit(-0.5 + 1.0 * u) else 0
        y = 1 if rng.random() < BASE * RR_A ** a * RR_UY ** u else 0
        yield {"id": i, "exposed": a, "event": y}


def bias_parameters():
    e1, e0 = expit(0.5), expit(-0.5)  # P(A = 1 | U = 1), P(A = 1 | U = 0)
    p1 = P_U * e1 / (P_U * e1 + (1 - P_U) * e0)
    p0 = P_U * (1 - e1) / (P_U * (1 - e1) + (1 - P_U) * (1 - e0))
    return {"p1": p1, "p0": p0, "rr_uy": RR_UY}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "exposed", "event"])
        writer.writeheader()
        writer.writerows(data)
    b = bias_parameters()
    print(f"wrote fixture.csv: {len(data)} people; causal log RR {math.log(RR_A):.10f}; "
          f"p1 {b['p1']:.6f}, p0 {b['p0']:.6f}, RR_UY {RR_UY}")
