# Fisher's exact test with the conditional maximum-likelihood odds ratio and its exact 95% interval
#
# The Python equivalent of r.R. scipy.stats.fisher_exact's statistic is the SAMPLE odds ratio ad/bc,
# not the conditional MLE R reports, so the odds ratio and its interval come from
# scipy.stats.contingency.odds_ratio(kind="conditional"); fisher_exact supplies only the p-value.

import math

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats.contingency import odds_ratio

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any()

a = int(((d.exposed == 1) & (d.outcome == 1)).sum()); b = int(((d.exposed == 1) & (d.outcome == 0)).sum())
c = int(((d.exposed == 0) & (d.outcome == 1)).sum()); dd = int(((d.exposed == 0) & (d.outcome == 0)).sum())

# PINNED: exposed in row 1, outcome in column 1, and kind="conditional" (the default kind, named).
table = np.array([[a, b], [c, dd]])
res = odds_ratio(table, kind="conditional")
ci = res.confidence_interval(confidence_level=0.95)
p = stats.fisher_exact(table, alternative="two-sided").pvalue

print(table)
print(f"conditional MLE odds ratio {res.statistic:.3f}, exact 95% {ci.low:.3f} to {ci.high:.3f}, p = {p:.4g}")
print(f"sample odds ratio ad/bc {a * dd / (b * c):.3f} (fisher_exact's statistic)")

print("\n--- HARNESS ---")
print(f"odds_ratio={res.statistic:.10f}")
print(f"log_odds_ratio_lcl={math.log(ci.low):.10f}")
print(f"log_odds_ratio_ucl={math.log(ci.high):.10f}")
print(f"log_odds_ratio={math.log(res.statistic):.10f}")
print(f"sample_odds_ratio={a * dd / (b * c):.10f}")
print(f"p_value={p:.10e}")
print(f"n={len(d)}")
