"""Decision curve analysis: the net benefit of acting on a prediction model's risk.

The same analysis as r.R: net benefit TP/n - FP/n x t / (1 - t) for treating when the predicted risk
exceeds t, beside treat-everyone (prevalence - (1 - prevalence) x t / (1 - t)) and treat-no-one (0), at
the prespecified thresholds 0.1, 0.2 and 0.3. See r.R for how to read it.
"""

import pandas as pd

d = pd.read_csv("fixture.csv")
n = len(d)
prev = d["event"].mean()


def nb(t):
    treat = d["predicted_risk"] > t
    return ((treat & (d.event == 1)).sum() - (treat & (d.event == 0)).sum() * t / (1 - t)) / n


def nb_all(t):
    return prev - (1 - prev) * t / (1 - t)


for t in (0.1, 0.2, 0.3):
    print(f"threshold {t:.1f}: model {nb(t):.4f}, treat all {nb_all(t):.4f}, treat none 0")

print("\n--- HARNESS ---")
for t in (0.1, 0.2, 0.3):
    key = f"{t:.1f}".replace(".", "")
    print(f"net_benefit_model_t{key}={nb(t):.10f}\nnet_benefit_all_t{key}={nb_all(t):.10f}")
print(f"prevalence={prev:.10f}\nn={n}")
