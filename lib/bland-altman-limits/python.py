"""Method comparison: Bland-Altman bias and 95% limits of agreement, with their intervals.

The same analysis as r.R: the bias (mean difference, B - A) with a t interval, the limits bias +/- 1.96 SD
with Bland and Altman's approximate intervals (SE sqrt(3 s^2 / n), t on n - 1 df), the proportional-bias
slope of differences on means, and the correlation, which is reported to show it is not agreement.
"""

import numpy as np
import pandas as pd
from scipy import stats

d = pd.read_csv("fixture.csv")
dif = (d["method_b"] - d["method_a"]).to_numpy()
avg = ((d["method_a"] + d["method_b"]) / 2).to_numpy()
n = len(dif)
bias, s = dif.mean(), dif.std(ddof=1)
tq, zq = stats.t.ppf(0.975, n - 1), stats.norm.ppf(0.975)
lo, hi = bias - zq * s, bias + zq * s
se_bias, se_loa = s / np.sqrt(n), np.sqrt(3 * s ** 2 / n)
slope = np.polyfit(avg, dif, 1)[0]
r = np.corrcoef(d["method_a"], d["method_b"])[0, 1]

print(f"bias {bias:.2f} (95% CI {bias - tq * se_bias:.2f} to {bias + tq * se_bias:.2f}); "
      f"limits of agreement {lo:.2f} to {hi:.2f}; correlation {r:.3f}")

print("\n--- HARNESS ---")
print(f"bias={bias:.10f}\nbias_lcl={bias - tq * se_bias:.10f}\nbias_ucl={bias + tq * se_bias:.10f}")
print(f"sd_difference={s:.10f}")
print(f"loa_lower={lo:.10f}\nloa_lower_lcl={lo - tq * se_loa:.10f}\nloa_lower_ucl={lo + tq * se_loa:.10f}")
print(f"loa_upper={hi:.10f}\nloa_upper_lcl={hi - tq * se_loa:.10f}\nloa_upper_ucl={hi + tq * se_loa:.10f}")
print(f"proportional_bias_slope={slope:.10f}\ncorrelation={r:.10f}\nn={n}")
