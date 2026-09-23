"""Statistical process control: u-, p- and individuals charts with baseline limits.

The same charts as r.R (R uses qcc), written out: centre lines from the baseline months 1 to 12; 3-sigma
limits applied to all 24 months, varying with each month's denominator for the u- and p-charts (the
p-chart's lower limit floored at 0); the individuals chart's sigma as the average moving range over the
baseline divided by d2 = 1.128. A signal is a point beyond a limit. See r.R for choosing the chart.
"""

import numpy as np
import pandas as pd

d = pd.read_csv("fixture.csv")
base = (d.month <= 12).to_numpy()

# u-chart: infections per catheter-day.
u_stat = (d.infections / d.catheter_days).to_numpy()
u_centre = d.infections[base].sum() / d.catheter_days[base].sum()
u_se = np.sqrt(u_centre / d.catheter_days.to_numpy())
u_ucl, u_lcl = u_centre + 3 * u_se, np.maximum(u_centre - 3 * u_se, 0)
u_beyond = np.flatnonzero((u_stat > u_ucl) | (u_stat < u_lcl)) + 1

# p-chart: readmissions per discharge.
p_stat = (d.readmissions / d.discharges).to_numpy()
p_centre = d.readmissions[base].sum() / d.discharges[base].sum()
p_se = np.sqrt(p_centre * (1 - p_centre) / d.discharges.to_numpy())
p_ucl, p_lcl = np.minimum(p_centre + 3 * p_se, 1), np.maximum(p_centre - 3 * p_se, 0)
p_beyond = np.flatnonzero((p_stat > p_ucl) | (p_stat < p_lcl)) + 1

# Individuals chart: mean length of stay.
x = d.mean_los.to_numpy()
i_centre = x[base].mean()
i_sigma = np.abs(np.diff(x[base])).mean() / 1.128
i_ucl, i_lcl = i_centre + 3 * i_sigma, i_centre - 3 * i_sigma
i_beyond = np.flatnonzero((x > i_ucl) | (x < i_lcl)) + 1

print(f"u-chart: centre {1000 * u_centre:.3f} per 1,000 catheter-days, signals in months {', '.join(map(str, u_beyond))}")
print(f"p-chart: centre {p_centre:.4f}, signals {len(p_beyond)}; I-chart: centre {i_centre:.3f}, signals {len(i_beyond)}")

print("\n--- HARNESS ---")
print(f"u_centre_per_1000={1000 * u_centre:.10f}\nu_ucl_month24_per_1000={1000 * u_ucl[-1]:.10f}")
print(f"u_signals={len(u_beyond)}\nu_first_signal_month={u_beyond[0] if len(u_beyond) else 0}")
print(f"p_centre={p_centre:.10f}\np_ucl_month24={p_ucl[-1]:.10f}\np_lcl_month24={p_lcl[-1]:.10f}\np_signals={len(p_beyond)}")
print(f"i_centre={i_centre:.10f}\ni_sigma={i_sigma:.10f}\ni_ucl={i_ucl:.10f}\ni_lcl={i_lcl:.10f}\ni_signals={len(i_beyond)}")
