"""Noninferiority analysis of a risk difference against a prespecified margin.

The same analysis as r.R: the risk difference (new minus standard) with Newcombe's hybrid score interval
(statsmodels' confint_proportions_2indep, method="newcomb"), the Wald lower limit beside it, and the
noninferiority conclusion (lower limit above the margin of -0.10) in the intention-to-treat and the
per-protocol populations. See r.R for why both are needed.
"""

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.proportion import confint_proportions_2indep

d = pd.read_csv("fixture.csv")
margin = -0.10
z = stats.norm.ppf(0.975)


def analyse(s):
    new, std = s[s.randomized_new == 1], s[s.randomized_new == 0]
    x1, n1, x0, n0 = int(new.success.sum()), len(new), int(std.success.sum()), len(std)
    p1, p0 = x1 / n1, x0 / n0
    lcl, ucl = confint_proportions_2indep(x1, n1, x0, n0, method="newcomb", compare="diff")
    wald = p1 - p0 - z * np.sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
    return p1 - p0, lcl, ucl, wald


itt = analyse(d)
pp = analyse(d[d.per_protocol == 1])
for name, r in (("ITT", itt), ("PP", pp)):
    print(f"{name}: difference {r[0]:.3f} (95% CI {r[1]:.3f} to {r[2]:.3f}); noninferior at margin {margin:.2f}: {r[1] > margin}")

print("\n--- HARNESS ---")
print(f"risk_difference_itt={itt[0]:.10f}\nitt_lcl={itt[1]:.10f}\nitt_ucl={itt[2]:.10f}\nitt_wald_lcl={itt[3]:.10f}")
print(f"risk_difference_pp={pp[0]:.10f}\npp_lcl={pp[1]:.10f}\npp_ucl={pp[2]:.10f}")
print(f"noninferior_itt={int(itt[1] > margin)}\nnoninferior_pp={int(pp[1] > margin)}")
print(f"n={len(d)}\nn_per_protocol={int(d.per_protocol.sum())}")
