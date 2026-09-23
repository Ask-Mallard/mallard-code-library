"""Seeded fixture for propensity-score matching: a confounded cohort with an effect that varies by severity.

Severity L ~ N(0, 1) and a binary comorbidity Z ~ Bernoulli(0.4) both raise the chance of treatment,
logit P(A = 1) = -0.8 + 0.6 L + 0.6 Z, and both raise the outcome. The treatment effect grows with
severity: y = 10 + (2 + 2.0 L) A + 1.5 L + 1.0 Z + noise (SD 2).

Because the effect varies with L and the treated are sicker, the three usual targets differ:
the ATE (everyone) is exactly 2, the ATT (the treated) is 2 + 2.0 E[L | A = 1], larger. Matching each
treated patient to a control estimates the ATT, so the fixture's truth is the ATT, integrated from the
data-generating process by `truth_att()`. A caliper that drops a few treated patients moves the target
slightly towards the ATE; the recovery tolerance is calibrated with that included.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261101
N = 2000
PS_INTERCEPT, PS_L, PS_Z = -0.8, 0.6, 0.6
EFFECT, EFFECT_BY_L = 2.0, 2.0
P_Z = 0.4


def expit(x):
    return 1.0 / (1.0 + math.exp(-x))


def propensity(l, z):
    return expit(PS_INTERCEPT + PS_L * l + PS_Z * z)


def rows(seed=SEED):
    rng = random.Random(seed)
    for i in range(1, N + 1):
        l = rng.gauss(0.0, 1.0)
        z = 1 if rng.random() < P_Z else 0
        a = 1 if rng.random() < propensity(l, z) else 0
        y = 10 + (EFFECT + EFFECT_BY_L * l) * a + 1.5 * l + 1.0 * z + rng.gauss(0.0, 2.0)
        yield {"id": i, "severity": round(l, 6), "comorbid": z, "treated": a, "y": round(y, 6)}


def _integrate(f, lo=-10.0, hi=10.0, steps=20000):
    # Simpson's rule over the standard normal density; the integrands are smooth and decay fast.
    h = (hi - lo) / steps
    total = 0.0
    for k in range(steps + 1):
        x = lo + k * h
        w = 1 if k in (0, steps) else (4 if k % 2 else 2)
        total += w * f(x) * math.exp(-0.5 * x * x) / math.sqrt(2 * math.pi)
    return total * h / 3


def truth_att():
    """E[effect | A = 1] = 2 + 2.0 E[L | A = 1], averaging over Z and L."""
    num = den = 0.0
    for z, pz in ((0, 1 - P_Z), (1, P_Z)):
        num += pz * _integrate(lambda l: propensity(l, z) * (EFFECT + EFFECT_BY_L * l))
        den += pz * _integrate(lambda l: propensity(l, z))
    return num / den


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "severity", "comorbid", "treated", "y"])
        writer.writeheader()
        writer.writerows(data)
    treated = sum(r["treated"] for r in data)
    print(f"wrote fixture.csv: {len(data)} patients, {treated} treated; ATT {truth_att():.10f}, ATE {EFFECT}")
