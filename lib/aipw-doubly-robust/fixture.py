"""Seeded fixture for augmented inverse-probability weighting (AIPW): a confounded cohort, binary outcome.

Severity L ~ N(0, 1) and a comorbidity Z ~ Bernoulli(0.4) raise both the chance of treatment,
logit P(A = 1) = -0.5 + 0.8 L + 0.5 Z, and the risk of the outcome,
logit P(Y = 1) = -1.0 + 0.7 A + 0.8 L + 0.5 Z.

The target is the average treatment effect as a marginal risk difference, E[Y(1)] - E[Y(0)], which
`truth()` integrates from the data-generating process. Both working models in the entry (logistic
propensity on L and Z, logistic outcome on A, L and Z) are correctly specified here; the control file
breaks each in turn to show that AIPW needs only one of them right.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261104
N = 2500
P_Z = 0.4
PS = (-0.5, 0.8, 0.5)
OUT = (-1.0, 0.7, 0.8, 0.5)  # intercept, treatment, L, Z


def expit(x):
    return 1.0 / (1.0 + math.exp(-x))


def risk(a, l, z):
    return expit(OUT[0] + OUT[1] * a + OUT[2] * l + OUT[3] * z)


def rows(seed=SEED):
    rng = random.Random(seed)
    for i in range(1, N + 1):
        l = rng.gauss(0.0, 1.0)
        z = 1 if rng.random() < P_Z else 0
        a = 1 if rng.random() < expit(PS[0] + PS[1] * l + PS[2] * z) else 0
        y = 1 if rng.random() < risk(a, l, z) else 0
        yield {"id": i, "severity": round(l, 6), "comorbid": z, "treated": a, "event": y}


def _integrate(f, lo=-10.0, hi=10.0, steps=20000):
    h = (hi - lo) / steps
    total = 0.0
    for k in range(steps + 1):
        x = lo + k * h
        w = 1 if k in (0, steps) else (4 if k % 2 else 2)
        total += w * f(x) * math.exp(-0.5 * x * x) / math.sqrt(2 * math.pi)
    return total * h / 3


def truth():
    """Marginal risks under treatment for everyone and for no one, and their contrasts."""
    r1 = r0 = 0.0
    for z, pz in ((0, 1 - P_Z), (1, P_Z)):
        r1 += pz * _integrate(lambda l: risk(1, l, z))
        r0 += pz * _integrate(lambda l: risk(0, l, z))
    return {"risk_treated": r1, "risk_untreated": r0, "risk_difference": r1 - r0,
            "log_risk_ratio": math.log(r1 / r0),
            "log_marginal_odds_ratio": math.log(r1 / (1 - r1)) - math.log(r0 / (1 - r0))}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "severity", "comorbid", "treated", "event"])
        writer.writeheader()
        writer.writerows(data)
    t = truth()
    print(f"wrote fixture.csv: {len(data)} patients, {sum(r['treated'] for r in data)} treated, "
          f"{sum(r['event'] for r in data)} events; RD {t['risk_difference']:.10f}, "
          f"log RR {t['log_risk_ratio']:.10f}, log marginal OR {t['log_marginal_odds_ratio']:.10f} "
          f"(conditional {OUT[1]})")
