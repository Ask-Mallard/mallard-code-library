"""Seeded fixture for a 1:2 matched case-control study.

WHY THE CSV IS COMMITTED AND NOT REGENERATED PER LANGUAGE. The obvious design is to seed each
language's own RNG and let it build the data. That would be wrong, and quietly so: R's Mersenne
Twister stream and NumPy's are different, so "seed 42" produces two different
datasets. Two implementations analysing two datasets cannot be compared to each other at all, and
the cross-language check -- the whole reason this library exists -- would be measuring sampling
variation while appearing to measure agreement.

So this script runs once, writes fixture.csv, and the CSV is committed. Every language reads
byte-identical rows. The generator is committed too, so the data is reproducible rather than
mysterious, but nothing in CI regenerates it.

WHY THE TRUE VALUE IS KNOWN EXACTLY. Within a matched set, choosing which member becomes the case
with probability proportional to exp(b1*x1 + b2*x2) IS the conditional likelihood that conditional
logistic regression maximises. The betas below are therefore the parameters of the model being
fitted, not merely of some data-generating process that resembles it, and clogit should recover
them up to sampling error. That is what makes a recovery tolerance meaningful rather than decorative.

Stdlib only, on purpose: this must run anywhere, including an environment with no scientific stack.
"""

import csv
import math
import random
from pathlib import Path

SEED = 20260906
N_SETS = 400          # matched sets
CONTROLS_PER_CASE = 2

# The true parameters. exposure log-OR of 0.6931 is an odds ratio of exactly 2.
BETA_EXPOSURE = math.log(2.0)
BETA_COVARIATE = 0.4


def main() -> None:
    rng = random.Random(SEED)
    rows = []

    for s in range(1, N_SETS + 1):
        members = []
        # A stratum-level baseline. Matching removes it from the conditional likelihood entirely,
        # which is the point of the design, so its value must not affect the recovered betas.
        # Leaving it large is a deliberate test of that: if an implementation silently fits an
        # UNconditional model, this term biases it and the recovery check fails.
        stratum_effect = rng.gauss(0.0, 2.0)

        for _ in range(1 + CONTROLS_PER_CASE):
            exposed = 1 if rng.random() < 0.35 else 0
            covariate = rng.gauss(0.0, 1.0)
            linpred = BETA_EXPOSURE * exposed + BETA_COVARIATE * covariate
            members.append({"exposed": exposed, "covariate": covariate, "linpred": linpred,
                            "stratum_effect": stratum_effect})

        # Choose the case proportional to exp(linpred): the 1:m conditional likelihood.
        weights = [math.exp(m["linpred"]) for m in members]
        total = sum(weights)
        draw = rng.random() * total
        case_index, running = 0, 0.0
        for i, w in enumerate(weights):
            running += w
            if draw <= running:
                case_index = i
                break

        for i, m in enumerate(members):
            rows.append({
                "set_id": s,
                "case": 1 if i == case_index else 0,
                "exposed": m["exposed"],
                "covariate": round(m["covariate"], 6),
            })

    out = Path(__file__).with_name("fixture.csv")
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["set_id", "case", "exposed", "covariate"])
        writer.writeheader()
        writer.writerows(rows)

    cases = sum(r["case"] for r in rows)
    exposed_cases = sum(r["exposed"] for r in rows if r["case"] == 1)
    exposed_controls = sum(r["exposed"] for r in rows if r["case"] == 0)
    print(f"wrote {out.name}: {len(rows)} rows, {N_SETS} sets, {cases} cases")
    print(f"  exposed among cases:    {exposed_cases}/{cases}")
    print(f"  exposed among controls: {exposed_controls}/{len(rows) - cases}")
    print(f"  true exposure OR: {math.exp(BETA_EXPOSURE):.4f}   true covariate beta: {BETA_COVARIATE}")


if __name__ == "__main__":
    main()
