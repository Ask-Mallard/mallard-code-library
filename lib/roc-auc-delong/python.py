"""ROC analysis: the area under the curve with a DeLong interval, and a threshold.

The same analysis as r.R, written out: the AUC as the Mann-Whitney probability (ties count one half),
DeLong's variance from the placement values of each case and control, a normal interval on the AUC
scale (as pROC's ci.auc with method = "delong"), and the threshold maximizing Youden's index among
midpoints between consecutive observed values, as pROC's coords. Direction pinned: higher values mean
disease. See r.R for what each part does and does not say.
"""

import numpy as np
import pandas as pd
from scipy import stats

d = pd.read_csv("fixture.csv")
x = d.loc[d.disease == 1, "marker"].to_numpy()
y = d.loc[d.disease == 0, "marker"].to_numpy()
m, n = len(x), len(y)

# PINNED: higher marker = disease. psi(case, control) = 1, 1/2 or 0.
psi = (x[:, None] > y[None, :]) + 0.5 * (x[:, None] == y[None, :])
auc = psi.mean()
v10 = psi.mean(axis=1)  # placement of each case
v01 = psi.mean(axis=0)  # placement of each control
se = float(np.sqrt(v10.var(ddof=1) / m + v01.var(ddof=1) / n))
z = stats.norm.ppf(0.975)

values = np.unique(d["marker"].to_numpy())
cuts = np.r_[-np.inf, (values[:-1] + values[1:]) / 2, np.inf]
sens = np.array([(x > c).mean() for c in cuts])
spec = np.array([(y <= c).mean() for c in cuts])
k = int(np.argmax(sens + spec))

print(f"AUC {auc:.3f} (DeLong 95% CI {auc - z * se:.3f} to {auc + z * se:.3f}); Youden threshold {cuts[k]:.3f} "
      f"(sensitivity {sens[k]:.3f}, specificity {spec[k]:.3f})")

print("\n--- HARNESS ---")
print(f"auc={auc:.10f}\nauc_se={se:.10f}\nauc_lcl={auc - z * se:.10f}\nauc_ucl={auc + z * se:.10f}")
print(f"youden_threshold={cuts[k]:.10f}\nsensitivity={sens[k]:.10f}\nspecificity={spec[k]:.10f}")
print(f"cases={m}\ncontrols={n}")
