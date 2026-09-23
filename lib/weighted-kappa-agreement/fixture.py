"""Seeded fixture for agreement between two raters on an ordinal scale: Cohen's and weighted kappa.

300 patients are graded 1 to 4 by two radiologists. Each patient's true grade T has probabilities
0.30, 0.30, 0.25, 0.15. Each rater, independently given T, reports T with probability 0.70 and otherwise
a neighbouring grade (split equally between the two neighbours; at the ends, the one neighbour). Most
disagreements are therefore one grade apart, which quadratic weighting treats as near-agreement.

`truth()` computes the population kappa (unweighted) and quadratic-weighted kappa exactly from the joint
distribution of the two ratings.

Stdlib only. `rows(seed)` is the generator; running the file writes the committed CSV.
"""

import csv
import random

SEED = 20261123
N = 300
P_TRUE = (0.30, 0.30, 0.25, 0.15)
P_EXACT = 0.70
K = 4


def p_report(r, t):
    """P(rater reports r | true grade t)."""
    if r == t:
        return P_EXACT
    neighbours = [g for g in (t - 1, t + 1) if 1 <= g <= K]
    return (1 - P_EXACT) / len(neighbours) if r in neighbours else 0.0


def _report(rng, t):
    u = rng.random()
    for r in range(1, K + 1):
        u -= p_report(r, t)
        if u < 0:
            return r
    return K


def rows(seed=SEED):
    rng = random.Random(seed)
    for i in range(1, N + 1):
        u, t = rng.random(), 1
        while u >= P_TRUE[t - 1]:
            u -= P_TRUE[t - 1]
            t += 1
        yield {"id": i, "rater_a": _report(rng, t), "rater_b": _report(rng, t)}


def truth():
    joint = [[sum(P_TRUE[t - 1] * p_report(a, t) * p_report(b, t) for t in range(1, K + 1))
              for b in range(1, K + 1)] for a in range(1, K + 1)]
    row = [sum(joint[a]) for a in range(K)]
    col = [sum(joint[a][b] for a in range(K)) for b in range(K)]
    out = {}
    for name, w in (("kappa", lambda a, b: 1.0 if a == b else 0.0),
                    ("kappa_quadratic", lambda a, b: 1 - (a - b) ** 2 / (K - 1) ** 2)):
        po = sum(w(a, b) * joint[a][b] for a in range(K) for b in range(K))
        pe = sum(w(a, b) * row[a] * col[b] for a in range(K) for b in range(K))
        out[name] = (po - pe) / (1 - pe)
    return out


if __name__ == "__main__":
    data = list(rows())
    with open("fixture.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "rater_a", "rater_b"])
        writer.writeheader()
        writer.writerows(data)
    t = truth()
    print(f"wrote fixture.csv: {len(data)} patients; kappa {t['kappa']:.10f}, quadratic {t['kappa_quadratic']:.10f}")
