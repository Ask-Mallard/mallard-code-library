# Kaplan-Meier estimation and the log-rank test
#
# The Python equivalent of the R file beside this one. lifelines rather than statsmodels: it is the
# maintained survival library and the one whose confidence bands can be compared with R's.

import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

d = pd.read_csv("fixture.csv")
unexposed = d[d["exposed"] == 0]
exposed = d[d["exposed"] == 1]

# PINNED DEFAULT 1: the exponential Greenwood interval, which IS the log-log transform R is told to
# use above. lifelines does this and offers no alternative, which happens to match the pin; that is
# worth stating rather than relying on, because R's DEFAULT is the plain log transform, so the same
# curve carries a different default band in each language.
kmf0, kmf1 = KaplanMeierFitter(), KaplanMeierFitter()
kmf0.fit(unexposed["time"], unexposed["event"], label="unexposed")
kmf1.fit(exposed["time"], exposed["event"], label="exposed")

# PINNED DEFAULT 2: the unweighted log-rank. lifelines' logrank_test is that test; the weighted
# family lives in other functions.
lr = logrank_test(exposed["time"], unexposed["time"],
                  exposed["event"], unexposed["event"])

print("Kaplan-Meier, by arm")
print(f"  unexposed median {kmf0.median_survival_time_}, exposed median {kmf1.median_survival_time_}")
print("\nlog-rank test")
lr.print_summary()

# `.asof` takes the last row at or before 12, which is what a step function's value AT 12 means.
# Indexing for an exact 12.0 would depend on an event having landed on that grid point.
ci0 = kmf0.confidence_interval_survival_function_.asof(12.0)

print("\n--- HARNESS ---")
print(f"median_unexposed={kmf0.median_survival_time_:.10f}")
print(f"median_exposed={kmf1.median_survival_time_:.10f}")
print(f"surv_unexposed_12={float(kmf0.predict(12.0)):.10f}")
print(f"surv_exposed_12={float(kmf1.predict(12.0)):.10f}")
print(f"ci_lower_unexposed_12={float(ci0.iloc[0]):.10f}")
print(f"ci_upper_unexposed_12={float(ci0.iloc[1]):.10f}")
print(f"logrank_chisq={lr.test_statistic:.10f}")
print(f"n_events={int(d['event'].sum())}")
print(f"n={len(d)}")
