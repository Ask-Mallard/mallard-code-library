# Table 1: baseline characteristics overall and by exposure group, with standardized differences
#
# The Python equivalent of r.R. The traps are defaults: numpy's std divides by n (ddof=0) while R's
# sd and pandas' std divide by n - 1, and quantile definitions differ between packages. Both are
# pinned below rather than trusted.

import numpy as np
import pandas as pd

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any() and set(d.exposed).issubset({0, 1})

e1 = d[d.exposed == 1]
e0 = d[d.exposed == 0]


def sd(x):
    # PINNED: ddof=1, the n - 1 divisor R's sd() uses. np.std(x) defaults to ddof=0.
    return float(np.std(x, ddof=1))


def smd_continuous(x1, x0):
    return float((x1.mean() - x0.mean()) / np.sqrt((sd(x1) ** 2 + sd(x0) ** 2) / 2))


def smd_binary(x1, x0):
    p1, p0 = x1.mean(), x0.mean()
    return float((p1 - p0) / np.sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2))


# PINNED: interpolation="linear", which is R's type 7.
q1, median, q3 = d.los_days.quantile([0.25, 0.5, 0.75], interpolation="linear")


def prop(x, level):
    return float((x == level).mean())


print(f"N = {len(d)} (exposed {len(e1)}, unexposed {len(e0)})")
print(f"Age, mean (SD)          {d.age.mean():.1f} ({sd(d.age):.1f}) | "
      f"{e1.age.mean():.1f} ({sd(e1.age):.1f}) vs {e0.age.mean():.1f} ({sd(e0.age):.1f}), "
      f"SMD {smd_continuous(e1.age, e0.age):.2f}")
print(f"Length of stay, median [IQR] {median:.2f} [{q1:.2f}, {q3:.2f}]")
print(f"Diabetes, n (%)        {int(d.diabetes.sum())} ({100 * d.diabetes.mean():.1f}) | "
      f"SMD {smd_binary(e1.diabetes, e0.diabetes):.2f}")
for level in ("never", "former", "current"):
    print(f"Smoking {level:<8} n (%)  {int((d.smoking == level).sum())} "
          f"({100 * prop(d.smoking, level):.1f})")

print("\n--- HARNESS ---")
out = {
    "n": len(d), "n_exposed": len(e1), "n_unexposed": len(e0),
    "age_mean": d.age.mean(), "age_sd": sd(d.age),
    "age_mean_exposed": e1.age.mean(), "age_mean_unexposed": e0.age.mean(),
    "age_sd_exposed": sd(e1.age), "age_sd_unexposed": sd(e0.age),
    "age_smd": smd_continuous(e1.age, e0.age),
    "los_q1": q1, "los_median": median, "los_q3": q3,
    "diabetes_prop": d.diabetes.mean(),
    "diabetes_prop_exposed": e1.diabetes.mean(), "diabetes_prop_unexposed": e0.diabetes.mean(),
    "diabetes_smd": smd_binary(e1.diabetes, e0.diabetes),
    "smoking_never_prop": prop(d.smoking, "never"),
    "smoking_former_prop": prop(d.smoking, "former"),
    "smoking_current_prop": prop(d.smoking, "current"),
}
for key, value in out.items():
    print(f"{key}={value}" if isinstance(value, int) else f"{key}={float(value):.10f}")
