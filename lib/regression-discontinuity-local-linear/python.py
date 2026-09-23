"""Sharp regression discontinuity: local linear estimation with robust bias-corrected inference.

The same analysis as r.R with the Python port of rdrobust, by the same authors (Calonico, Cattaneo,
Farrell and Titiunik). Every option is named rather than left to either port's defaults: local linear
(p = 1) with a local quadratic bias correction (q = 2), triangular kernel, one common MSE-optimal
bandwidth, nearest-neighbour variance with 3 matches. Report the conventional estimate with the ROBUST
bias-corrected interval; see r.R for why the conventional interval is too narrow at the MSE-optimal
bandwidth.
"""

import pandas as pd
import statsmodels.formula.api as smf
from rdrobust import rdrobust

d = pd.read_csv("fixture.csv")
r = rdrobust(y=d["y"], x=d["score"], c=0, p=1, q=2, kernel="triangular", bwselect="mserd",
             vce="nn", nnmatch=3)

coef, se, ci, bws = r.coef, r.se, r.ci, r.bws
est = float(coef.loc["Conventional"].iloc[0])
lcl, ucl = float(ci.loc["Robust"].iloc[0]), float(ci.loc["Robust"].iloc[1])

global_jump = smf.ols("y ~ treated * score", d).fit().params["treated"]

print(f"jump at the cut-off {est:.3f}; robust 95% CI {lcl:.3f} to {ucl:.3f}; bandwidth {bws.loc['h', 'left']:.3f}")

print("\n--- HARNESS ---")
print(f"jump={est:.10f}\njump_se_conventional={float(se.loc['Conventional'].iloc[0]):.10f}")
print(f"jump_bias_corrected={float(coef.loc['Bias-Corrected'].iloc[0]):.10f}\njump_se_robust={float(se.loc['Robust'].iloc[0]):.10f}")
print(f"robust_lcl={lcl:.10f}\nrobust_ucl={ucl:.10f}")
print(f"conventional_lcl={float(ci.loc['Conventional'].iloc[0]):.10f}\nconventional_ucl={float(ci.loc['Conventional'].iloc[1]):.10f}")
print(f"bandwidth_h={float(bws.loc['h', 'left']):.10f}\nbandwidth_b={float(bws.loc['b', 'left']):.10f}")
print(f"n_left_in_h={int(r.N_h[0])}\nn_right_in_h={int(r.N_h[1])}")
print(f"global_linear_jump={global_jump:.10f}\nn={len(d)}")
