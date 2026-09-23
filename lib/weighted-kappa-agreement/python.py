"""Agreement between two raters on an ordinal scale: Cohen's kappa and quadratic-weighted kappa.

The same analysis as r.R, written out: kappa = (observed - chance agreement) / (1 - chance), unweighted
and with quadratic weights 1 - (i - j)^2 / (k - 1)^2, each with the large-sample standard error of
Fleiss, Cohen and Everitt (1969). See r.R for why an ordinal scale needs weights and why percent
agreement is not enough.
"""

import numpy as np
import pandas as pd
from scipy import stats

d = pd.read_csv("fixture.csv")
lev = np.arange(1, 5)
k, n = len(lev), len(d)
p = np.zeros((k, k))
for a, b in zip(d["rater_a"], d["rater_b"]):
    p[a - 1, b - 1] += 1 / n
pr, pc = p.sum(axis=1), p.sum(axis=0)


def kappa_se(w):
    po, pe = np.sum(w * p), np.sum(w * np.outer(pr, pc))
    kap = (po - pe) / (1 - pe)
    wr, wc = w @ pc, w.T @ pr
    inner = np.sum(p * (w - (wr[:, None] + wc[None, :]) * (1 - kap)) ** 2)
    return kap, float(np.sqrt((inner - (kap - pe * (1 - kap)) ** 2) / (n * (1 - pe) ** 2)))


k0, se0 = kappa_se(np.eye(k))
k2, se2 = kappa_se(1 - (lev[:, None] - lev[None, :]) ** 2 / (k - 1) ** 2)

z = stats.norm.ppf(0.975)
print(f"percent agreement {np.trace(p):.3f}; kappa {k0:.3f} ({k0 - z * se0:.3f} to {k0 + z * se0:.3f}); "
      f"quadratic-weighted kappa {k2:.3f} ({k2 - z * se2:.3f} to {k2 + z * se2:.3f})")

print("\n--- HARNESS ---")
print(f"percent_agreement={np.trace(p):.10f}")
print(f"kappa={k0:.10f}\nkappa_se={se0:.10f}")
print(f"kappa_quadratic={k2:.10f}\nkappa_quadratic_se={se2:.10f}")
print(f"n={n}")
