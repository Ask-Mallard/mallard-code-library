"""Cochran-Armitage test for trend in a proportion across ordered groups, with the linear slope.

No standard Python package gives the Cochran-Armitage statistic in the form R's prop.trend.test
reports, so it is written out. statsmodels' linear-by-linear association test is a DIFFERENT
statistic, (N - 1) r^2 rather than N r^2; harness/categorical_examples.test.py shows they differ.
"""

import numpy as np
import pandas as pd
from scipy import stats

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any() and set(d.outcome).issubset({0, 1})

g = d.groupby("dose").outcome.agg(["sum", "count"]).sort_index()
score = g.index.to_numpy(dtype=float)
events, totals = g["sum"].to_numpy(dtype=float), g["count"].to_numpy(dtype=float)
big_n = totals.sum()
p_bar = events.sum() / big_n

# PINNED: scores are the dose values, and the statistic divides by N (Cochran-Armitage, as in R's
# prop.trend.test), not N - 1.
s_bar = (totals * score).sum() / big_n
numerator = (score * (events - totals * p_bar)).sum() ** 2
denominator = p_bar * (1 - p_bar) * (totals * (score - s_bar) ** 2).sum()
chi2 = numerator / denominator
p = stats.chi2.sf(chi2, 1)

# Weighted least-squares slope of the group proportions on the scores, weights = group sizes.
prop = events / totals
slope = (totals * (score - s_bar) * (prop - p_bar)).sum() / (totals * (score - s_bar) ** 2).sum()

print(g.assign(proportion=prop))
print(f"trend chi-square {chi2:.4f} on 1 df, p = {p:.4g}; slope {slope:.4f} per step")

print("\n--- HARNESS ---")
print(f"trend_chi2={chi2:.10f}")
print(f"p_value={p:.10e}")
print(f"slope={slope:.10f}")
print(f"n={int(big_n)}")
