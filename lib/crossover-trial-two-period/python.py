# Two-period crossover (AB/BA): the period-adjusted treatment effect
#
# The Python equivalent of r.R, using scipy's two-sample t-test on the within-patient period
# differences, independent of R's t.test.

import pandas as pd
from scipy import stats

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any()

w = d.pivot(index=["id", "sequence"], columns="period", values="y").reset_index()
w["diff"] = w[1] - w[2]

# PINNED: period differences compared BY SEQUENCE (BA first), halved; equal_var=True.
ba, ab = w["diff"][w.sequence == "BA"], w["diff"][w.sequence == "AB"]
fit = stats.ttest_ind(ba, ab, equal_var=True)
ci = fit.confidence_interval(confidence_level=0.95)
effect = (ba.mean() - ab.mean()) / 2

print(f"treatment effect (B - A) {effect:.3f}, 95% {ci.low / 2:.3f} to {ci.high / 2:.3f}, {int(fit.df)} df")

print("\n--- HARNESS ---")
print(f"treatment_effect={effect:.10f}\ntreatment_effect_lcl={ci.low / 2:.10f}\ntreatment_effect_ucl={ci.high / 2:.10f}")
print(f"df={int(fit.df)}\npatients={len(w)}")
