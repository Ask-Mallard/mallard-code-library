"""Seeded fixture for rate regression on overdispersed counts with unequal follow-up.

Counts are drawn from a GAMMA-POISSON mixture: each patient's expected count is
`person_time * exp(Xb)` multiplied by a gamma frailty with mean 1 and variance 0.5. That is exactly
the negative binomial (NB2) model, so the coefficients below are the parameters BOTH a Poisson and
a negative binomial regression estimate -- the Poisson coefficients stay consistent under
overdispersion, and only their standard errors go wrong.

TWO PROPERTIES ARE DELIBERATE.

FOLLOW-UP VARIES, from six months to five years. With equal follow-up the offset can be omitted
entirely and the rate ratio comes out right anyway, so a fixture with constant person-time cannot
test the offset -- and the offset is the thing most easily got wrong here, because an offset must be
the LOG of the person-time and passing the person-time itself fits a different model that still runs.

THE COUNTS ARE OVERDISPERSED, with a gamma variance of 0.5. Under a true Poisson the naive and the
robust standard errors converge, so an equidispersed fixture could not tell them apart and the
whole point of the entry would be untestable. With this frailty the naive Poisson standard error is
materially too small, which is the failure that matters: the estimate is right and the interval is
too narrow, so the analysis is confidently wrong rather than visibly wrong.

Stdlib only. The Poisson sampler is Knuth's, which is exact and fine at these means.
"""

import csv
import math
import random
from pathlib import Path

SEED = 20260908
N = 3000
LOG_BASE_RATE = math.log(0.30)      # events per person-year in the unexposed at covariate 0
BETA_EXPOSURE = math.log(1.6)       # a true RATE RATIO of exactly 1.6
BETA_COVARIATE = 0.3
GAMMA_VARIANCE = 0.5                # NB2 dispersion: theta = 1 / 0.5 = 2
FOLLOWUP_MIN = 0.5
FOLLOWUP_MAX = 5.0


def poisson(rng: random.Random, mu: float) -> int:
    """Knuth's algorithm. Exact, and the means here are small enough that it never underflows."""
    limit, k, p = math.exp(-mu), 0, 1.0
    while True:
        p *= rng.random()
        if p <= limit:
            return k
        k += 1


def main() -> None:
    rng = random.Random(SEED)
    rows = []
    for i in range(1, N + 1):
        exposed = 1 if rng.random() < 0.45 else 0
        covariate = rng.gauss(0.0, 1.0)
        person_time = rng.uniform(FOLLOWUP_MIN, FOLLOWUP_MAX)
        mu = person_time * math.exp(LOG_BASE_RATE + BETA_EXPOSURE * exposed
                                    + BETA_COVARIATE * covariate)
        # Gamma frailty with mean 1 and variance GAMMA_VARIANCE. gammavariate takes (shape, scale).
        frailty = rng.gammavariate(1.0 / GAMMA_VARIANCE, GAMMA_VARIANCE)
        rows.append({"id": i, "events": poisson(rng, mu * frailty),
                     "person_time": round(person_time, 6), "exposed": exposed,
                     "covariate": round(covariate, 6)})

    out = Path(__file__).with_name("fixture.csv")
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["id", "events", "person_time", "exposed", "covariate"])
        w.writeheader()
        w.writerows(rows)

    total = sum(r["events"] for r in rows)
    pt = sum(r["person_time"] for r in rows)
    print(f"wrote {out.name}: {len(rows)} rows, {total} events over {pt:.1f} person-years")
    print(f"  follow-up ranges {min(r['person_time'] for r in rows):.3f} to "
          f"{max(r['person_time'] for r in rows):.3f} years -- THE OFFSET MATTERS")

    # Two crude checks needing no model.
    rates = {}
    for arm in (0, 1):
        sub = [r for r in rows if r["exposed"] == arm]
        ev = sum(r["events"] for r in sub)
        t = sum(r["person_time"] for r in sub)
        rates[arm] = ev / t
        print(f"  arm {arm}: {len(sub):4d} patients, {ev:5d} events, {t:8.1f} person-years, "
              f"rate {ev / t:.4f}/year")
    print(f"  crude rate ratio {rates[1] / rates[0]:.4f}   true "
          f"{math.exp(BETA_EXPOSURE):.4f}")

    mean = total / len(rows)
    var = sum((r["events"] - mean) ** 2 for r in rows) / (len(rows) - 1)
    print(f"  counts: mean {mean:.4f}, variance {var:.4f}, ratio {var / mean:.4f}")
    print("  A RATIO ABOVE 1 IS THE OVERDISPERSION. Under a true Poisson it would be about 1 and")
    print("  the naive and robust standard errors would agree, testing nothing.")


if __name__ == "__main__":
    main()
