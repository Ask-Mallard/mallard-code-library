"""Predictive values and likelihood ratios from a diagnostic accuracy study.

The same analysis as r.R: PPV and NPV as single proportions with Wilson score intervals (statsmodels'
method="wilson", which is prop.test without continuity correction), likelihood ratios with log-scale
intervals (Simel, Samsa and Matchar 1991), and the PPV recomputed at a prevalence of 0.05. See r.R for
why predictive values belong to a population and likelihood ratios to the test.
"""

import math

import pandas as pd
from scipy import stats
from statsmodels.stats.proportion import proportion_confint

d = pd.read_csv("fixture.csv")
tp = int(((d.disease == 1) & (d.test_positive == 1)).sum())
fn = int(((d.disease == 1) & (d.test_positive == 0)).sum())
fp = int(((d.disease == 0) & (d.test_positive == 1)).sum())
tn = int(((d.disease == 0) & (d.test_positive == 0)).sum())
nd, nnd = tp + fn, fp + tn

ppv, ppv_ci = tp / (tp + fp), proportion_confint(tp, tp + fp, method="wilson")
npv, npv_ci = tn / (tn + fn), proportion_confint(tn, tn + fn, method="wilson")
se, sp = tp / nd, tn / nnd
lr_pos, lr_neg = se / (1 - sp), (1 - se) / sp
se_log_pos = math.sqrt((1 - se) / (se * nd) + sp / ((1 - sp) * nnd))
se_log_neg = math.sqrt(se / ((1 - se) * nd) + (1 - sp) / (sp * nnd))
p_other = 0.05
ppv_other = se * p_other / (se * p_other + (1 - sp) * (1 - p_other))

z = stats.norm.ppf(0.975)
print(f"PPV {ppv:.3f} ({ppv_ci[0]:.3f} to {ppv_ci[1]:.3f}), NPV {npv:.3f}; LR+ {lr_pos:.2f}, LR- {lr_neg:.3f}; "
      f"PPV at 5% prevalence {ppv_other:.3f}")

print("\n--- HARNESS ---")
print(f"ppv={ppv:.10f}\nppv_lcl={ppv_ci[0]:.10f}\nppv_ucl={ppv_ci[1]:.10f}")
print(f"npv={npv:.10f}\nnpv_lcl={npv_ci[0]:.10f}\nnpv_ucl={npv_ci[1]:.10f}")
print(f"log_lr_positive={math.log(lr_pos):.10f}\nlog_lr_positive_se={se_log_pos:.10f}")
print(f"log_lr_negative={math.log(lr_neg):.10f}\nlog_lr_negative_se={se_log_neg:.10f}")
print(f"ppv_at_prevalence_0.05={ppv_other:.10f}")
print(f"n={len(d)}")
