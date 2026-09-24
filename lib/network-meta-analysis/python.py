"""Network meta-analysis (frequentist, random effects) of three treatments from two-arm trials.

The same model as r.R (netmeta), written out for two-arm trials, where netmeta's graph-theoretical
estimator is weighted least squares on the trial contrasts: each contrast treat1 - treat2 is a row of
the design matrix over the basic parameters (B and C against the reference A). The common-effect fit
gives Q; tau^2 is the generalised DerSimonian-Laird estimate (Q - df) / (tr(W) - tr(W X (X'WX)^-1 X'W)),
floored at 0; the random-effects fit reweights by 1 / (se^2 + tau^2). The between-design (inconsistency)
component of Q is Q minus the within-design Q from fitting each design separately. See r.R for the
assumptions a network needs.
"""

import numpy as np
import pandas as pd
from scipy import stats

d = pd.read_csv("fixture.csv")
params = ["B", "C"]  # basic parameters against the reference A
X = np.zeros((len(d), 2))
for i, (t1, t2) in enumerate(zip(d.treat1, d.treat2)):
    if t1 in params:
        X[i, params.index(t1)] += 1
    if t2 in params:
        X[i, params.index(t2)] -= 1
y, s2 = d.TE.to_numpy(float), d.seTE.to_numpy(float) ** 2


def wls(w):
    info = (X * w[:, None]).T @ X
    beta = np.linalg.solve(info, (X * w[:, None]).T @ y)
    return beta, np.linalg.inv(info)


w = 1 / s2
beta_c, _ = wls(w)
q = float(np.sum(w * (y - X @ beta_c) ** 2))
df = len(y) - X.shape[1]
H = (X * w[:, None]) @ np.linalg.inv((X * w[:, None]).T @ X) @ (X * w[:, None]).T
tau2 = max(0.0, (q - df) / (np.sum(w) - np.trace(H)))
beta, cov = wls(1 / (s2 + tau2))

# Within-design heterogeneity: each design's own common-effect Q; the rest of Q is between designs.
designs = d.treat1 + ":" + d.treat2
q_within = 0.0
for g in designs.unique():
    m = (designs == g).to_numpy()
    mu = np.sum(w[m] * y[m]) / np.sum(w[m])
    q_within += np.sum(w[m] * (y[m] - mu) ** 2)
q_incons = q - q_within

contrast = {"B_vs_A": np.array([1.0, 0.0]), "C_vs_A": np.array([0.0, 1.0]), "C_vs_B": np.array([-1.0, 1.0])}
est = {k: float(c @ beta) for k, c in contrast.items()}
se = {k: float(np.sqrt(c @ cov @ c)) for k, c in contrast.items()}
z = stats.norm.ppf(0.975)

for k in ("B_vs_A", "C_vs_A"):
    print(f"{k.replace('_vs_', ' vs ')}: OR {np.exp(est[k]):.3f} (95% CI {np.exp(est[k] - z * se[k]):.3f} to "
          f"{np.exp(est[k] + z * se[k]):.3f})")
print(f"C vs B: OR {np.exp(est['C_vs_B']):.3f}; tau^2 {tau2:.4f}; Q {q:.3f} (between designs {q_incons:.3f})")

print("\n--- HARNESS ---")
for k in ("B_vs_A", "C_vs_A", "C_vs_B"):
    print(f"log_or_{k}={est[k]:.10f}\nlog_or_{k}_se={se[k]:.10f}")
print(f"common_log_or_C_vs_B={float(np.array([-1.0, 1.0]) @ beta_c):.10f}")
print(f"tau2={tau2:.10f}\nq={q:.10f}\nq_inconsistency={q_incons:.10f}")
print(f"studies={len(d)}")
