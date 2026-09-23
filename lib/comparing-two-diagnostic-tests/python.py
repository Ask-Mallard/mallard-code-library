"""Comparing two diagnostic tests on the same patients: paired AUCs and paired sensitivities.

The same analysis as r.R, written out. The AUCs and their DeLong covariance come from each marker's
placement values over the same cases and controls; the paired standard error of the difference uses the
covariance, and the unpaired one (ignoring it) is reported beside it. At the prespecified thresholds, the
sensitivities are compared among patients with disease by McNemar's test (no continuity correction) with
a Wald interval for the paired difference. Direction pinned: higher values mean disease.
"""

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.contingency_tables import mcnemar

d = pd.read_csv("fixture.csv")
cases, controls = d[d.disease == 1], d[d.disease == 0]
m, n0 = len(cases), len(controls)


def placements(col):
    x, y = cases[col].to_numpy(), controls[col].to_numpy()
    psi = (x[:, None] > y[None, :]) + 0.5 * (x[:, None] == y[None, :])
    return psi.mean(), psi.mean(axis=1), psi.mean(axis=0)


auc_a, v10a, v01a = placements("marker_a")
auc_b, v10b, v01b = placements("marker_b")
S10 = np.cov(np.vstack([v10a, v10b]))
S01 = np.cov(np.vstack([v01a, v01b]))
S = S10 / m + S01 / n0
diff = auc_a - auc_b
se_paired = float(np.sqrt(S[0, 0] + S[1, 1] - 2 * S[0, 1]))
se_unpaired = float(np.sqrt(S[0, 0] + S[1, 1]))

pa = (cases["marker_a"] > 0.6).astype(int).to_numpy()
pb = (cases["marker_b"] > 0.6).astype(int).to_numpy()
b01, b10 = int(((pa == 1) & (pb == 0)).sum()), int(((pa == 0) & (pb == 1)).sum())
table = [[int(((pa == 0) & (pb == 0)).sum()), b10], [b01, int(((pa == 1) & (pb == 1)).sum())]]
mc = mcnemar(table, exact=False, correction=False)
sens_diff = (b01 - b10) / m
sens_diff_se = np.sqrt((b01 + b10) - (b01 - b10) ** 2 / m) / m

zstat = diff / se_paired
p = 2 * stats.norm.sf(abs(zstat))
print(f"AUC A {auc_a:.3f}, B {auc_b:.3f}, difference {diff:.3f} (paired DeLong p = {p:.4f}); "
      f"sensitivity difference {sens_diff:.3f} (McNemar p = {mc.pvalue:.4f})")

print("\n--- HARNESS ---")
print(f"auc_a={auc_a:.10f}\nauc_b={auc_b:.10f}\nauc_difference={diff:.10f}")
print(f"auc_difference_se={se_paired:.10f}\nauc_difference_se_unpaired={se_unpaired:.10f}\nauc_test_z={zstat:.10f}")
print(f"sensitivity_a={pa.mean():.10f}\nsensitivity_b={pb.mean():.10f}\nsensitivity_difference={sens_diff:.10f}\n"
      f"sensitivity_difference_se={sens_diff_se:.10f}")
print(f"mcnemar_statistic={mc.statistic:.10f}\ndiscordant_a_only={b01}\ndiscordant_b_only={b10}\nn={len(d)}")
