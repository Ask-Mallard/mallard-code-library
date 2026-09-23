"""Seeded fixture for causal mediation with a binary mediator and a binary outcome.

A covariate C ~ N(0, 1) confounds everything. Treatment A, a binary mediator M and a binary outcome Y:

logit P(A = 1) = -0.3 + 0.5 C
logit P(M = 1) = -0.5 + 1.0 A + 0.4 C
logit P(Y = 1) = -1.5 + 0.4 A + 1.2 M + 0.3 A M + 0.5 C

The natural direct and indirect effects on the risk-difference scale are defined by counterfactuals,
E[Y(a, M(a*))], and `truth()` integrates them from the data-generating process. Because the outcome
model is logistic, the product-of-coefficients and difference-of-coefficients methods do not estimate
them: the odds ratio is non-collapsible, and those methods have no term for the A x M interaction.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261110
N = 4000
G = (-0.5, 1.0, 0.4)             # mediator model: intercept, A, C
B = (-1.5, 0.4, 1.2, 0.3, 0.5)   # outcome model: intercept, A, M, A x M, C


def expit(x):
    return 1.0 / (1.0 + math.exp(-x))


def p_m(a, c):
    return expit(G[0] + G[1] * a + G[2] * c)


def p_y(a, m, c):
    return expit(B[0] + B[1] * a + B[2] * m + B[3] * a * m + B[4] * c)


def rows(seed=SEED):
    rng = random.Random(seed)
    for i in range(1, N + 1):
        c = rng.gauss(0.0, 1.0)
        a = 1 if rng.random() < expit(-0.3 + 0.5 * c) else 0
        m = 1 if rng.random() < p_m(a, c) else 0
        y = 1 if rng.random() < p_y(a, m, c) else 0
        yield {"id": i, "c": round(c, 6), "treated": a, "mediator": m, "event": y}


def _integrate(f, lo=-10.0, hi=10.0, steps=20000):
    h = (hi - lo) / steps
    total = 0.0
    for k in range(steps + 1):
        x = lo + k * h
        w = 1 if k in (0, steps) else (4 if k % 2 else 2)
        total += w * f(x) * math.exp(-0.5 * x * x) / math.sqrt(2 * math.pi)
    return total * h / 3


def potential_risk(a, a_star):
    """E[Y(a, M(a*))]: the outcome risk under treatment a with the mediator as it would be under a*."""
    return _integrate(lambda c: p_y(a, 1, c) * p_m(a_star, c) + p_y(a, 0, c) * (1 - p_m(a_star, c)))


def truth():
    r11, r10, r00 = potential_risk(1, 1), potential_risk(1, 0), potential_risk(0, 0)
    return {"natural_direct_effect": r10 - r00, "natural_indirect_effect": r11 - r10, "total_effect": r11 - r00}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "c", "treated", "mediator", "event"])
        writer.writeheader()
        writer.writerows(data)
    t = truth()
    print(f"wrote fixture.csv: {len(data)} participants; NDE {t['natural_direct_effect']:.10f}, "
          f"NIE {t['natural_indirect_effect']:.10f}, total {t['total_effect']:.10f}")
