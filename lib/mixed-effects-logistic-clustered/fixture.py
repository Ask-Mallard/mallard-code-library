"""Seeded fixture for a cluster-randomised binary outcome with a random intercept.

Exposure is assigned AT THE CLUSTER LEVEL, which is the case where clustering actually bites and
the cluster-randomized trial with a binary outcome. Individual-level
exposure inside clusters would be a much weaker test: the cluster random effect would then be
nearly orthogonal to it and an analysis ignoring clustering would still land near the truth.

The outcome is drawn from logit(p) = b0 + b1*exposed + u_j with u_j ~ N(0, sigma^2). That is a
random-intercept logistic model by construction, so b1 is the CLUSTER-CONDITIONAL log odds ratio
that glmer estimates.

THE CONDITIONAL PARAMETER IS NOT THE MARGINAL ONE. A population-average model (GEE) targets a
different, attenuated quantity on the same data, and for a nonlinear link the two are genuinely
different numbers rather than different spellings. sigma is set high enough that the gap is
visible, because an entry whose fixture made the two coincide could not demonstrate the distinction
it exists to teach.

Stdlib only.
"""

import csv
import math
import random
from pathlib import Path

SEED = 20260906
N_CLUSTERS = 40
PER_CLUSTER = 30
BETA0 = -1.0
BETA_EXPOSURE = math.log(2.0)   # true CLUSTER-CONDITIONAL odds ratio of exactly 2
SIGMA_U = 0.8                   # between-cluster SD on the logit scale


def main() -> None:
    rng = random.Random(SEED)
    rows = []
    for c in range(1, N_CLUSTERS + 1):
        exposed = 1 if c % 2 == 0 else 0        # half the clusters, allocated at cluster level
        u = rng.gauss(0.0, SIGMA_U)
        for k in range(PER_CLUSTER):
            eta = BETA0 + BETA_EXPOSURE * exposed + u
            p = 1.0 / (1.0 + math.exp(-eta))
            rows.append({"cluster": c, "exposed": exposed,
                         "outcome": 1 if rng.random() < p else 0})

    out = Path(__file__).with_name("fixture.csv")
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["cluster", "exposed", "outcome"])
        w.writeheader()
        w.writerows(rows)

    icc = SIGMA_U ** 2 / (SIGMA_U ** 2 + math.pi ** 2 / 3)
    n1 = [r for r in rows if r["exposed"] == 1]
    n0 = [r for r in rows if r["exposed"] == 0]
    p1 = sum(r["outcome"] for r in n1) / len(n1)
    p0 = sum(r["outcome"] for r in n0) / len(n0)
    crude_or = (p1 / (1 - p1)) / (p0 / (1 - p0))
    print(f"wrote {out.name}: {len(rows)} rows, {N_CLUSTERS} clusters of {PER_CLUSTER}")
    print(f"  latent-scale ICC: {icc:.3f}   between-cluster SD: {SIGMA_U}")
    print(f"  crude MARGINAL odds ratio: {crude_or:.4f}")
    print(f"  true CLUSTER-CONDITIONAL odds ratio: {math.exp(BETA_EXPOSURE):.4f}")
    print("  the gap between those two is the point: they are different estimands, and the")
    print("  marginal one is attenuated toward 1 under a nonlinear link.")


if __name__ == "__main__":
    main()
