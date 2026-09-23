# McNemar's test for a paired binary outcome, with the change in proportion and its 95% interval
#
# The Python equivalent of r.R, using statsmodels' mcnemar (independent of R's mcnemar.test and
# binom.test) for both tests. The change and its paired Wald interval are a short formula written in
# both files.

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.contingency_tables import mcnemar

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any()

n = len(d)
b = int(((d.before == 1) & (d.after == 0)).sum())
c = int(((d.before == 0) & (d.after == 1)).sum())
table = np.array([[int(((d.before == 1) & (d.after == 1)).sum()), b],
                  [c, int(((d.before == 0) & (d.after == 0)).sum())]])

# PINNED: exact=False, correction=False for the chi-square (statsmodels' default is the EXACT test),
# and exact=True for the binomial p-value.
chi = mcnemar(table, exact=False, correction=False)
exact = mcnemar(table, exact=True)

change = (c - b) / n
se = np.sqrt(b + c - (c - b) ** 2 / n) / n
z = stats.norm.ppf(0.975)

print(f"discordant: {b} positive->negative, {c} negative->positive")
print(f"change in proportion positive {change:.4f}, 95% {change - z * se:.4f} to {change + z * se:.4f}")
print(f"McNemar chi-square {chi.statistic:.3f}; exact p {exact.pvalue:.4g}")

print("\n--- HARNESS ---")
print(f"change={change:.10f}\nchange_lcl={change - z * se:.10f}\nchange_ucl={change + z * se:.10f}")
print(f"discordant_pos_neg={b}\ndiscordant_neg_pos={c}")
print(f"mcnemar_chi2={chi.statistic:.10f}")
print(f"exact_p_value={exact.pvalue:.10e}")
print(f"n={n}")
