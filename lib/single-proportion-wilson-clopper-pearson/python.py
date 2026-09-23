# A single proportion (a prevalence) with Wilson and Clopper-Pearson 95% intervals
#
# The Python equivalent of r.R. statsmodels' proportion_confint is an implementation independent of
# R's binom package and of base R's binom.test, so the two languages agreeing is evidence rather
# than the same arithmetic typed twice.

import pandas as pd
from statsmodels.stats.proportion import proportion_confint

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any() and set(d.outcome).issubset({0, 1})

x = int(d.outcome.sum())
n = len(d)

# PINNED DEFAULT: method="wilson". proportion_confint defaults to "normal", the Wald interval,
# which is too narrow near 0 and 1 and can run below 0.
wilson_lcl, wilson_ucl = proportion_confint(x, n, alpha=0.05, method="wilson")
# "beta" is statsmodels' name for Clopper-Pearson. There is no method called "exact".
cp_lcl, cp_ucl = proportion_confint(x, n, alpha=0.05, method="beta")

print(f"events {x} of {n}, proportion {x / n:.4f}")
print(f"Wilson 95%          {wilson_lcl:.4f} to {wilson_ucl:.4f}")
print(f"Clopper-Pearson 95% {cp_lcl:.4f} to {cp_ucl:.4f}")

print("\n--- HARNESS ---")
print(f"prevalence={x / n:.10f}")
print(f"wilson_lcl={wilson_lcl:.10f}")
print(f"wilson_ucl={wilson_ucl:.10f}")
print(f"cp_lcl={cp_lcl:.10f}")
print(f"cp_ucl={cp_ucl:.10f}")
print(f"events={x}")
print(f"n={n}")
