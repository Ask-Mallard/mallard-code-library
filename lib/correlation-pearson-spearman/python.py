# Pearson and Spearman correlation, each with a 95% interval
#
# The Python equivalent of r.R. scipy's pearsonr gives the Fisher z interval directly (independent of
# R's cor.test); spearmanr gives no interval, so it is built with the same Fieller standard error.

import numpy as np
import pandas as pd
from scipy import stats

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any()
n = len(d)
z = stats.norm.ppf(0.975)

pear = stats.pearsonr(d.bmi, d.sbp)
pci = pear.confidence_interval(confidence_level=0.95)

rho = stats.spearmanr(d.bmi, d.sbp).statistic
# PINNED: the Fieller standard error 1.06 / sqrt(n - 3) for Spearman.
se_s = 1.06 / np.sqrt(n - 3)
sp_lcl, sp_ucl = np.tanh(np.arctanh(rho) + np.array([-1, 1]) * z * se_s)

print(f"Pearson r {pear.statistic:.4f}, 95% {pci.low:.4f} to {pci.high:.4f}")
print(f"Spearman rho {rho:.4f}, 95% {sp_lcl:.4f} to {sp_ucl:.4f}")

print("\n--- HARNESS ---")
print(f"pearson_r={pear.statistic:.10f}\npearson_lcl={pci.low:.10f}\npearson_ucl={pci.high:.10f}")
print(f"spearman_rho={rho:.10f}\nspearman_lcl={sp_lcl:.10f}\nspearman_ucl={sp_ucl:.10f}")
print(f"n={n}")
