"""Random-effects meta-analysis of odds ratios: REML for tau^2, Hartung-Knapp for the interval.

The same analysis as r.R (metafor), written out: log odds ratios with Woolf variances; tau^2 by
maximizing the restricted (REML) log-likelihood over tau^2 >= 0 to 1e-12; the mean with the
Hartung-Knapp interval (variance rescaled by the weighted residual sum of squares over k - 1, t on k - 1
df); I^2 from metafor's typical within-trial variance; the prediction interval on k - 1 df; and
DerSimonian-Laird with a z interval beside it. statsmodels' combine_effects offers DL and Paule-Mandel
but not REML with Hartung-Knapp. See r.R for what each quantity means.
"""

import numpy as np
import pandas as pd
from scipy import optimize, stats

d = pd.read_csv("fixture.csv")
a, c = d.events_treated.to_numpy(float), d.events_control.to_numpy(float)
b, dd = d.n_treated.to_numpy(float) - a, d.n_control.to_numpy(float) - c
y = np.log(a * dd / (b * c))
v = 1 / a + 1 / b + 1 / c + 1 / dd
k = len(y)


def pooled(tau2):
    w = 1 / (v + tau2)
    mu = np.sum(w * y) / np.sum(w)
    return w, mu


def neg_restricted_loglik(tau2):
    w, mu = pooled(tau2)
    return 0.5 * (np.sum(np.log(v + tau2)) + np.log(np.sum(w)) + np.sum(w * (y - mu) ** 2))


res = optimize.minimize_scalar(neg_restricted_loglik, bounds=(0, 10), method="bounded", options={"xatol": 1e-14})
tau2 = max(res.x, 0.0) if neg_restricted_loglik(res.x) < neg_restricted_loglik(0.0) else 0.0
w, mu = pooled(tau2)
k1 = k - 1
s2 = np.sum(w * (y - mu) ** 2) / k1          # Hartung-Knapp rescaling
se = np.sqrt(s2 / np.sum(w))
tq = stats.t.ppf(0.975, k1)
w0 = 1 / v
typical = k1 * np.sum(w0) / (np.sum(w0) ** 2 - np.sum(w0 ** 2))
i2 = 100 * tau2 / (tau2 + typical)
q = np.sum(w0 * (y - np.sum(w0 * y) / np.sum(w0)) ** 2)
pi_half = tq * np.sqrt(tau2 + se ** 2)

tau2_dl = max(0.0, (q - k1) / (np.sum(w0) - np.sum(w0 ** 2) / np.sum(w0)))
w_dl, mu_dl = pooled(tau2_dl)

print(f"summary OR {np.exp(mu):.3f} (HK 95% CI {np.exp(mu - tq * se):.3f} to {np.exp(mu + tq * se):.3f}); "
      f"tau^2 {tau2:.4f}, I^2 {i2:.1f}%; prediction interval {np.exp(mu - pi_half):.3f} to {np.exp(mu + pi_half):.3f}")

print("\n--- HARNESS ---")
print(f"log_or={mu:.10f}\nlog_or_se={se:.10f}\nlog_or_lcl={mu - tq * se:.10f}\nlog_or_ucl={mu + tq * se:.10f}")
print(f"tau2={tau2:.10f}\ni2={i2:.10f}\nq={q:.10f}")
print(f"prediction_lcl={mu - pi_half:.10f}\nprediction_ucl={mu + pi_half:.10f}")
print(f"dl_log_or={mu_dl:.10f}\ndl_se_z={np.sqrt(1 / np.sum(w_dl)):.10f}\ntrials={k}")
