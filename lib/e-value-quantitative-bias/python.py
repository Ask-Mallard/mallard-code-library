"""Sensitivity to unmeasured confounding: the E-value and a simple quantitative bias analysis.

The same analysis as r.R. Python has no E-value package, so the formula is written out
(VanderWeele and Ding 2017): for a risk ratio above 1, E = RR + sqrt(RR (RR - 1)), computed for the point
estimate and for the confidence limit closer to 1 (1 if the interval includes it). The bias-adjusted
risk ratio divides the observed one by (p1 (RR_UY - 1) + 1) / (p0 (RR_UY - 1) + 1). See r.R for what each
means and what the bias factor assumes.
"""

import math

import pandas as pd
from scipy import stats

d = pd.read_csv("fixture.csv")
a = int(d.loc[d.exposed == 1, "event"].sum()); n1 = int((d.exposed == 1).sum())
c0 = int(d.loc[d.exposed == 0, "event"].sum()); n0 = int((d.exposed == 0).sum())
log_rr = math.log((a / n1) / (c0 / n0))
se = math.sqrt(1 / a - 1 / n1 + 1 / c0 - 1 / n0)
z = stats.norm.ppf(0.975)
rr, lo, hi = math.exp(log_rr), math.exp(log_rr - z * se), math.exp(log_rr + z * se)


def e_value(r):
    r = r if r >= 1 else 1 / r  # a protective RR is inverted first
    return r + math.sqrt(r * (r - 1))


# The limit closer to the null: the lower one for RR > 1, the upper one for RR < 1; 1 if it spans 1.
near = lo if rr > 1 else hi
ev_point = e_value(rr)
ev_ci = 1.0 if (lo <= 1 <= hi) else e_value(near)

# PINNED: the bias parameters (external evidence in practice; the fixture's true values here).
p1, p0, rr_uy = 0.5236, 0.2879, 2.5
bias = (p1 * (rr_uy - 1) + 1) / (p0 * (rr_uy - 1) + 1)

print(f"observed RR {rr:.2f} (95% CI {lo:.2f} to {hi:.2f}); E-value {ev_point:.2f} (CI limit {ev_ci:.2f}); "
      f"bias-adjusted RR {rr / bias:.2f}")

print("\n--- HARNESS ---")
print(f"log_rr_observed={log_rr:.10f}\nlog_rr_observed_se={se:.10f}")
print(f"e_value_point={ev_point:.10f}\ne_value_ci={ev_ci:.10f}")
print(f"bias_factor={bias:.10f}\nlog_rr_bias_adjusted={log_rr - math.log(bias):.10f}")
print(f"n={len(d)}")
