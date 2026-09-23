"""Seeded fixture for target trial emulation with a grace period: clone, censor, weight.

An observational cohort of 3,000 eligible patients is followed for 10 intervals. In each interval k a
time-varying severity L_k is observed; an untreated patient may start treatment only in intervals 0, 1
or 2 (the grace period), with probability expit(-1 + 0.8 L_k), so sicker patients start sooner. Then an
event occurs with probability expit(-3 + 0.6 L_k - 0.7 A_k), where A_k = 1 once treatment has started.
Severity evolves as L_(k+1) = 0.8 L_k - 0.3 A_k + N(0, 0.5^2); L_0 ~ N(0, 1).

The target trial compares two strategies over the 10 intervals:
  (1) start treatment within the grace period: as observed in intervals 0 and 1, and in interval 2 for
      everyone who has not yet started;
  (2) never start.
`truth()` simulates 400,000 patients under each strategy (a fixed seed; Monte Carlo SD of the risk
difference about 0.001) and returns the 10-interval risks and their difference.

Two naive analyses are biased: comparing patients who ever started with those who never did gives
immortal time to the treated (they had to survive until they started), and comparing by the treatment
actually received ignores that sicker patients started.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import math
import random

SEED = 20261114
N = 3000
INTERVALS = 10
GRACE = 3
N_TRUTH = 400000
TRUTH_SEED = 99


def expit(x):
    return 1.0 / (1.0 + math.exp(-x))


def p_start(l):
    return expit(-1.0 + 0.8 * l)


def p_event(l, a):
    return expit(-3.0 + 0.6 * l - 0.7 * a)


def _follow(rng, strategy):
    """One patient's history. strategy: None (observed), 'grace' or 'never'. Returns a list of intervals."""
    l = rng.gauss(0.0, 1.0)
    a = 0
    out = []
    for k in range(INTERVALS):
        if not a and k < GRACE:
            start = rng.random() < p_start(l)  # drawn under every strategy, so the streams stay aligned
            if strategy == "never":
                start = False
            elif strategy == "grace" and k == GRACE - 1:
                start = True
            a = 1 if start else 0
        y = 1 if rng.random() < p_event(l, a) else 0
        out.append((k, round(l, 6), a, y))
        if y:
            break
        l = 0.8 * l - 0.3 * a + rng.gauss(0.0, 0.5)
    return out


def rows(seed=SEED):
    rng = random.Random(seed)
    for i in range(1, N + 1):
        for k, l, a, y in _follow(rng, None):
            yield {"id": i, "interval": k, "severity": l, "treated": a, "event": y}


def truth():
    risks = {}
    for strategy in ("grace", "never"):
        rng = random.Random(TRUTH_SEED)
        risks[strategy] = sum(_follow(rng, strategy)[-1][3] for _ in range(N_TRUTH)) / N_TRUTH
    return {"risk_grace": risks["grace"], "risk_never": risks["never"],
            "risk_difference": risks["grace"] - risks["never"]}


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "interval", "severity", "treated", "event"])
        writer.writeheader()
        writer.writerows(data)
    t = truth()
    print(f"wrote fixture.csv: {N} patients, {len(data)} patient-intervals, {sum(r['event'] for r in data)} events; "
          f"risk grace {t['risk_grace']:.6f}, never {t['risk_never']:.6f}, difference {t['risk_difference']:.10f}")
