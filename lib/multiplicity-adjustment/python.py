"""Multiplicity: adjusting ten p-values by Holm, Hochberg and Benjamini-Hochberg.

The same analysis as r.R: Welch's t-test on each of the ten outcomes, then statsmodels' multipletests
with method "holm", "simes-hochberg" (Hochberg's step-up) and "fdr_bh", which return the same adjusted
p-values as R's p.adjust. See r.R for which error rate each controls and when each is valid.
"""

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

d = pd.read_csv("fixture.csv")
t, c = d[d.treated == 1], d[d.treated == 0]
p = np.array([stats.ttest_ind(t[f"y{j}"], c[f"y{j}"], equal_var=False).pvalue for j in range(1, 11)])
diff1 = t["y1"].mean() - c["y1"].mean()
holm = multipletests(p, method="holm")[1]
hoch = multipletests(p, method="simes-hochberg")[1]
bh = multipletests(p, method="fdr_bh")[1]

print(f"rejected at 0.05: unadjusted {(p < 0.05).sum()}, Holm {(holm < 0.05).sum()}, "
      f"Hochberg {(hoch < 0.05).sum()}, Benjamini-Hochberg {(bh < 0.05).sum()}")

print("\n--- HARNESS ---")
for j in range(10):
    print(f"p_{j + 1}={p[j]:.10f}\np_holm_{j + 1}={holm[j]:.10f}\np_hochberg_{j + 1}={hoch[j]:.10f}\np_bh_{j + 1}={bh[j]:.10f}")
print(f"difference_1={diff1:.10f}")
print(f"rejected_unadjusted={(p < 0.05).sum()}\nrejected_holm={(holm < 0.05).sum()}\n"
      f"rejected_hochberg={(hoch < 0.05).sum()}\nrejected_bh={(bh < 0.05).sum()}")
