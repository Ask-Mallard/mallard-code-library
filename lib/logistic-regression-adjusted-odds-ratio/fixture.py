"""Seeded fixture for an ADJUSTED ODDS RATIO from a binary outcome.

The outcome is drawn from a LOGIT-link model, so the conditional log odds ratio the generator uses
is exactly the parameter logistic regression estimates. The truth is therefore a property of the
data rather than of any implementation of the analysis -- the same construction the conditional
logistic and modified Poisson entries use, applied to the logit link.

The covariate does two jobs, and both are deliberate:

  1. IT CONFOUNDS. Exposure depends on it, so the crude odds ratio is not the adjusted one and an
     entry about ADJUSTMENT has something to adjust for. An "adjusted" example whose covariate
     changed nothing would demonstrate the mechanics and none of the point.

  2. IT MAKES THE ODDS RATIO NON-COLLAPSIBLE, visibly. This generator computes the marginal odds
     ratio by standardising the TRUE model over the observed covariate distribution: no confounding
     remains in that comparison, and the marginal odds ratio is still not 2.0. That gap is not bias
     and not sampling error. It is what non-collapsibility means, in numbers, from the model that
     made the data.

Stdlib only, so it runs anywhere.
"""

import csv
import math
import random
from pathlib import Path

SEED = 20260915
N = 4000

BETA0 = -1.6
BETA_EXPOSURE = math.log(2.0)   # true CONDITIONAL odds ratio of exactly 2.0
BETA_COVARIATE = 0.5

# Exposure depends on the covariate, so the covariate confounds.
ALPHA0 = -0.3
ALPHA_COVARIATE = 0.8


def expit(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x))


def odds(p: float) -> float:
    return p / (1.0 - p)


def main() -> None:
    rng = random.Random(SEED)
    rows = []
    for i in range(1, N + 1):
        covariate = rng.gauss(0.0, 1.0)
        exposed = 1 if rng.random() < expit(ALPHA0 + ALPHA_COVARIATE * covariate) else 0
        p = expit(BETA0 + BETA_EXPOSURE * exposed + BETA_COVARIATE * covariate)
        rows.append({"id": i, "exposed": exposed, "covariate": round(covariate, 6),
                     "outcome": 1 if rng.random() < p else 0})

    out = Path(__file__).with_name("fixture.csv")
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["id", "exposed", "covariate", "outcome"])
        w.writeheader()
        w.writerows(rows)

    events = sum(r["outcome"] for r in rows)
    n1 = sum(1 for r in rows if r["exposed"] == 1)
    r1 = sum(r["outcome"] for r in rows if r["exposed"] == 1) / max(1, n1)
    r0 = sum(r["outcome"] for r in rows if r["exposed"] == 0) / max(1, N - n1)
    crude_or = odds(r1) / odds(r0)

    # THE MARGINAL ODDS RATIO, STANDARDISED OVER THE OBSERVED COVARIATE DISTRIBUTION, FROM THE TRUE
    # MODEL. Every row contributes both potential outcomes, so nothing here is confounded and
    # nothing is estimated: this is the generating model averaged, not a fit. It still does not
    # equal 2.0, and that is the non-collapsibility of the odds ratio rather than any error.
    cov = [r["covariate"] for r in rows]
    p1 = sum(expit(BETA0 + BETA_EXPOSURE + BETA_COVARIATE * z) for z in cov) / N
    p0 = sum(expit(BETA0 + BETA_COVARIATE * z) for z in cov) / N
    marginal_or = odds(p1) / odds(p0)
    marginal_rr = p1 / p0

    print(f"wrote {out.name}: {len(rows)} rows, {events} events ({100 * events / N:.1f}%)")
    print(f"  exposed: {n1} ({100 * n1 / N:.1f}%)")
    print(f"  CONDITIONAL odds ratio (what logistic regression estimates): {math.exp(BETA_EXPOSURE):.4f}")
    print(f"  crude odds ratio, unadjusted:                                {crude_or:.4f}   <- confounded")
    print(f"  MARGINAL odds ratio, standardised over the true model:       {marginal_or:.4f}   <- not 2.0, and not bias")
    print(f"  marginal risk ratio, same standardisation:                   {marginal_rr:.4f}   <- a different estimand again")


if __name__ == "__main__":
    main()
