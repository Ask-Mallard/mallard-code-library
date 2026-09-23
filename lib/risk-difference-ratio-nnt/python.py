# Risk difference, risk ratio and number needed to treat from a two-arm trial, with 95% intervals
#
# The Python equivalent of r.R, using statsmodels' confint_proportions_2indep for both intervals, an
# implementation independent of the Newcombe and Katz formulas written out in r.R.

import numpy as np
import pandas as pd
from statsmodels.stats.proportion import confint_proportions_2indep

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any()

x1 = int(d.outcome[d.treated == 1].sum()); n1 = int((d.treated == 1).sum())
x0 = int(d.outcome[d.treated == 0].sum()); n0 = int((d.treated == 0).sum())
p1, p0 = x1 / n1, x0 / n0

# PINNED: method="newcomb" for the difference (the default is "wald") and method="log" for the ratio
# (Katz). Treated is the FIRST count: compare="diff" returns p1 - p0.
rd_lcl, rd_ucl = confint_proportions_2indep(x1, n1, x0, n0, method="newcomb", compare="diff", alpha=0.05)
rr_lcl, rr_ucl = confint_proportions_2indep(x1, n1, x0, n0, method="log", compare="ratio", alpha=0.05)
rd = p1 - p0

# NNT bounds invert the RD bounds only when they exclude zero.
assert rd_ucl < 0 or rd_lcl > 0, "the RD interval crosses zero: report the NNT as two ranges"
nnt = 1 / abs(rd)
nnt_lcl, nnt_ucl = sorted(1 / abs(np.array([rd_lcl, rd_ucl])))

print(f"treated {x1}/{n1} ({p1:.3f}), control {x0}/{n0} ({p0:.3f})")
print(f"risk difference {rd:.4f}, Newcombe 95% {rd_lcl:.4f} to {rd_ucl:.4f}")
print(f"risk ratio {p1 / p0:.3f}, 95% {rr_lcl:.3f} to {rr_ucl:.3f}")
print(f"NNT {nnt:.1f}, 95% {nnt_lcl:.1f} to {nnt_ucl:.1f}")

print("\n--- HARNESS ---")
print(f"risk_treated={p1:.10f}\nrisk_control={p0:.10f}")
print(f"risk_difference={rd:.10f}\nrisk_difference_lcl={rd_lcl:.10f}\nrisk_difference_ucl={rd_ucl:.10f}")
print(f"log_risk_ratio={np.log(p1 / p0):.10f}")
print(f"risk_ratio={p1 / p0:.10f}\nrisk_ratio_lcl={rr_lcl:.10f}\nrisk_ratio_ucl={rr_ucl:.10f}")
print(f"nnt={nnt:.10f}\nnnt_lcl={nnt_lcl:.10f}\nnnt_ucl={nnt_ucl:.10f}")
print(f"n={len(d)}")
