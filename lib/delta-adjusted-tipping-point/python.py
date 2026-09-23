"""Delta-adjusted tipping-point analysis for missing outcomes (a pattern-mixture sensitivity analysis).

The same analysis as r.R: within each arm, the outcome is regressed on the baseline covariate among
patients with an outcome; each arm's mean averages observed and predicted outcomes, with a shift delta
added to the predictions for the TREATED patients whose outcome is missing; the standard error comes from
the stacked estimating equations. The tipping point is the delta at which the lower 95% limit reaches
zero, found by Brent's method as uniroot does. See r.R for what it means.
"""

import numpy as np
import pandas as pd
from scipy import optimize, stats

d = pd.read_csv("fixture.csv")
d["r"] = d["y"].notna().astype(int)


def arm_fit(a):
    s = d[d.treated == a].reset_index(drop=True)
    X = np.column_stack([np.ones(len(s)), s["x"]])
    obs = s["r"].to_numpy() == 1
    b = np.linalg.solve(X[obs].T @ X[obs], X[obs].T @ s["y"].to_numpy()[obs])
    return s, X, obs, b


arms = {1: arm_fit(1), 0: arm_fit(0)}


def effect_at(delta):
    est, var = {}, 0.0
    for a in (1, 0):
        s, X, obs, b = arms[a]
        pred = X @ b + (delta if a == 1 else 0.0)
        y = s["y"].to_numpy()
        val = np.where(obs, y, pred)
        mu = val.mean()
        res = np.where(obs, y - X @ b, 0.0)
        psi = np.column_stack([X * res[:, None], val - mu])
        A = np.zeros((3, 3))
        A[:2, :2] = -X[obs].T @ X[obs]
        A[2, :2] = (X * (~obs)[:, None]).sum(axis=0)
        A[2, 2] = -len(X)
        Ai = np.linalg.inv(A)
        V = Ai @ (psi.T @ psi) @ Ai.T
        est[a] = mu
        var += V[2, 2]
    return est[1] - est[0], np.sqrt(var)


z = stats.norm.ppf(0.975)
effect, se = effect_at(0.0)
tip = optimize.brentq(lambda dl: effect_at(dl)[0] - z * effect_at(dl)[1], -50, 0, xtol=1e-12)
cc = d.loc[d.treated == 1, "y"].mean() - d.loc[d.treated == 0, "y"].mean()

print(f"effect under MAR {effect:.3f} (95% CI {effect - z * se:.3f} to {effect + z * se:.3f}); "
      f"tipping point: treated missing outcomes {-tip:.2f} lower")

print("\n--- HARNESS ---")
print(f"effect={effect:.10f}\neffect_se={se:.10f}")
print(f"effect_lcl={effect - z * se:.10f}\neffect_ucl={effect + z * se:.10f}")
print(f"tipping_delta={tip:.10f}")
print(f"effect_at_delta_minus2={effect_at(-2.0)[0]:.10f}")
print(f"complete_case_difference={cc:.10f}\nmissing={int((d.r == 0).sum())}\nn={len(d)}")
