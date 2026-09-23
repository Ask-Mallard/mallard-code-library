"""Propensity-score matching: the average treatment effect in the TREATED.

The same analysis as r.R. Python has no matching package whose defaults are documented to match
MatchIt's, so the algorithm is written out, following MatchIt's nearest-neighbour matching step for step:

- distance: the LOGIT of a logistic propensity score (MatchIt's link = "linear.logit");
- caliper: 0.2 standard deviations of that distance, over the whole sample (std.caliper = TRUE);
- order: treated patients matched from the LARGEST distance down (m.order = "largest");
- 1:1, without replacement: each control is used at most once; a treated patient with no unused control
  inside the caliper is dropped;
- ties: the control that comes first in the data.

The effect is the matched difference in means; its SE clusters on the matched pair, as sandwich::vcovCL
does (statsmodels' cluster covariance applies the same G/(G-1) x (N-1)/(N-K) correction).
Matching estimates the ATT, not the ATE: see r.R.
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

d = pd.read_csv("fixture.csv")

ps = smf.glm("treated ~ severity + comorbid", d, family=sm.families.Binomial()).fit()
p = np.asarray(ps.fittedvalues)
dist = np.log(p / (1 - p))  # PINNED: the logit of the score, as link = "linear.logit"
caliper = 0.2 * np.std(dist, ddof=1)  # PINNED: 0.2 SD of the distance, whole sample

treated = np.flatnonzero(d["treated"].to_numpy() == 1)
controls = np.flatnonzero(d["treated"].to_numpy() == 0)
available = np.ones(len(controls), dtype=bool)
pairs = []
for t in treated[np.argsort(-dist[treated], kind="stable")]:  # PINNED: largest distance first
    gap = np.abs(dist[controls] - dist[t])
    gap[~available] = np.inf
    j = int(np.argmin(gap))  # first minimum: ties go to the control earliest in the data
    if gap[j] <= caliper:
        available[j] = False
        pairs.append((t, controls[j]))

rows = [(t, k) for k, (t, c) in enumerate(pairs)] + [(c, k) for k, (t, c) in enumerate(pairs)]
md = d.iloc[[r for r, _ in rows]].assign(subclass=[k for _, k in rows])

fit = smf.ols("y ~ treated", md).fit(cov_type="cluster", cov_kwds={"groups": md["subclass"]})
att, se = float(fit.params["treated"]), float(fit.bse["treated"])

sd_t = d.loc[d.treated == 1, "severity"].std(ddof=1)


def smd(x):
    return (x.loc[x.treated == 1, "severity"].mean() - x.loc[x.treated == 0, "severity"].mean()) / sd_t


naive = d.loc[d.treated == 1, "y"].mean() - d.loc[d.treated == 0, "y"].mean()
z = stats.norm.ppf(0.975)
print(f"ATT {att:.3f} (95% CI {att - z * se:.3f} to {att + z * se:.3f}); {len(pairs)} of {len(treated)} treated matched")

print("\n--- HARNESS ---")
print(f"att={att:.10f}\natt_se={se:.10f}\natt_lcl={att - z * se:.10f}\natt_ucl={att + z * se:.10f}")
print(f"naive_difference={naive:.10f}")
print(f"smd_severity_before={smd(d):.10f}\nsmd_severity_after={smd(md):.10f}")
print(f"matched_treated={len(pairs)}\nmatched_pairs={len(pairs)}\nn={len(d)}")
