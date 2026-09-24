"""Bootstrap confidence intervals for a ratio of means from skewed data: percentile and BCa.

The same procedure as r.R: a stratified bootstrap (patients resampled within each hospital), 2,000
resamples of the ratio of mean length of stay, the percentile interval and the BCa interval (bias
correction from the share of bootstrap values below the estimate, acceleration from the jackknife). The resamples come from the same shared
"minimal standard" generator as r.R, so both languages use identical resamples; in practice use
scipy.stats.bootstrap(method="BCa") with a seeded generator. Quantiles are numpy's default (R's type 7).
"""

import numpy as np
import pandas as pd
from scipy import stats

d = pd.read_csv("fixture.csv")
a = d.loc[d.hospital == "A", "los"].to_numpy()
b = d.loc[d.hospital == "B", "los"].to_numpy()
estimate = a.mean() / b.mean()
B = 2000

state = 20261201


def draw(n):
    """n uniform indices in 0..n-1 from the shared generator (Park-Miller, multiplier 48271)."""
    global state
    idx = np.empty(n, dtype=int)
    for k in range(n):
        state = (48271 * state) % 2147483647
        idx[k] = int(state / 2147483647 * n)
    return idx


boot_stat = np.array([a[draw(len(a))].mean() / b[draw(len(b))].mean() for _ in range(B)])
percentile = np.quantile(boot_stat, [0.025, 0.975])

z0 = stats.norm.ppf(np.mean(boot_stat < estimate))
jack = np.r_[[np.delete(a, i).mean() / b.mean() for i in range(len(a))],
             [a.mean() / np.delete(b, i).mean() for i in range(len(b))]]
jm = jack.mean()
acc = np.sum((jm - jack) ** 3) / (6 * np.sum((jm - jack) ** 2) ** 1.5)
zq = stats.norm.ppf([0.025, 0.975])
adj = stats.norm.cdf(z0 + (z0 + zq) / (1 - acc * (z0 + zq)))
bca = np.quantile(boot_stat, adj)

print(f"ratio of means {estimate:.3f}; percentile 95% CI {percentile[0]:.3f} to {percentile[1]:.3f}; "
      f"BCa {bca[0]:.3f} to {bca[1]:.3f}; bootstrap SE {boot_stat.std(ddof=1):.3f}")

print("\n--- HARNESS ---")
print(f"mean_ratio={estimate:.10f}\nbootstrap_se={boot_stat.std(ddof=1):.10f}")
print(f"percentile_lcl={percentile[0]:.10f}\npercentile_ucl={percentile[1]:.10f}")
print(f"bca_lcl={bca[0]:.10f}\nbca_ucl={bca[1]:.10f}\nbias_correction_z0={z0:.10f}\nacceleration={acc:.10f}")
print(f"resamples={B}")
