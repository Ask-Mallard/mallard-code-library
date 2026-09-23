"""Inter-rater reliability: the intraclass correlation, with the model, type and unit named.

The same analysis as r.R, written out from the two-way ANOVA mean squares (subjects, raters, residual):
ICC(A,1), two-way random, absolute agreement, single rater, with McGraw and Wong's F-based interval
(its Satterthwaite-type df), and ICC(C,1), consistency, with its F interval. The formulas follow
irr::icc line for line. See r.R for why the model, type and unit must be named.
"""

import numpy as np
import pandas as pd
from scipy import stats

d = pd.read_csv("fixture.csv")
Y = d.filter(like="rater_").to_numpy(float)
ns, nr = Y.shape
alpha = 0.05

ss_total = Y.var(ddof=1) * (ns * nr - 1)
ms_r = Y.mean(axis=1).var(ddof=1) * nr   # subjects
ms_c = Y.mean(axis=0).var(ddof=1) * ns   # raters
ms_e = (ss_total - ms_r * (ns - 1) - ms_c * (nr - 1)) / ((ns - 1) * (nr - 1))

# ICC(A,1) and its interval.
icc_a = (ms_r - ms_e) / (ms_r + (nr - 1) * ms_e + (nr / ns) * (ms_c - ms_e))
a = nr * icc_a / (ns * (1 - icc_a))
b = 1 + nr * icc_a * (ns - 1) / (ns * (1 - icc_a))
v = (a * ms_c + b * ms_e) ** 2 / ((a * ms_c) ** 2 / (nr - 1) + (b * ms_e) ** 2 / ((ns - 1) * (nr - 1)))
fl = stats.f.ppf(1 - alpha / 2, ns - 1, v)
fu = stats.f.ppf(1 - alpha / 2, v, ns - 1)
lo_a = ns * (ms_r - fl * ms_e) / (fl * (nr * ms_c + (nr * ns - nr - ns) * ms_e) + ns * ms_r)
hi_a = ns * (fu * ms_r - ms_e) / (nr * ms_c + (nr * ns - nr - ns) * ms_e + ns * fu * ms_r)

# ICC(C,1) and its interval.
icc_c = (ms_r - ms_e) / (ms_r + (nr - 1) * ms_e)
f_l = (ms_r / ms_e) / stats.f.ppf(1 - alpha / 2, ns - 1, (ns - 1) * (nr - 1))
f_u = (ms_r / ms_e) * stats.f.ppf(1 - alpha / 2, (ns - 1) * (nr - 1), ns - 1)
lo_c, hi_c = (f_l - 1) / (f_l + nr - 1), (f_u - 1) / (f_u + nr - 1)

print(f"ICC(A,1) {icc_a:.3f} (95% CI {lo_a:.3f} to {hi_a:.3f}); ICC(C,1) {icc_c:.3f}")

print("\n--- HARNESS ---")
print(f"icc_agreement={icc_a:.10f}\nicc_agreement_lcl={lo_a:.10f}\nicc_agreement_ucl={hi_a:.10f}")
print(f"icc_consistency={icc_c:.10f}\nicc_consistency_lcl={lo_c:.10f}\nicc_consistency_ucl={hi_c:.10f}")
print(f"subjects={ns}\nraters={nr}")
