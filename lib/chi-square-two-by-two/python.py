# Two-by-two table: chi-square, Fisher's exact test and the odds ratio
#
# The Python equivalent of the R file beside this one. Both languages apply Yates' correction by
# default and both are told not to.

import numpy as np
import pandas as pd
from scipy import stats

d = pd.read_csv("fixture.csv")

# Built in the same orientation as the R file: rows exposed then unexposed, columns outcome then no
# outcome. Written out cell by cell rather than with a crosstab, because a crosstab's row and column
# order follows the sort order of the values, which for 0/1 integers is the wrong way round.
a = int(((d["exposed"] == 1) & (d["outcome"] == 1)).sum())
b = int(((d["exposed"] == 1) & (d["outcome"] == 0)).sum())
c = int(((d["exposed"] == 0) & (d["outcome"] == 1)).sum())
e = int(((d["exposed"] == 0) & (d["outcome"] == 0)).sum())
tab = np.array([[a, b], [c, e]])

# PINNED DEFAULT: correction=False.
#
# scipy.stats.chi2_contingency applies Yates' continuity correction to any 2x2 by default, exactly
# as R's chisq.test does, so both files switch it off. On this fixture
# that is 49.90 uncorrected against 49.19 corrected, with every expected count in the hundreds.
pearson = stats.chi2_contingency(tab, correction=False)
yates = stats.chi2_contingency(tab, correction=True)
fisher = stats.fisher_exact(tab)

print("2x2 table (rows: exposed, unexposed; columns: outcome, no outcome)")
print(tab)
print(f"\nPearson chi-square, uncorrected: {pearson.statistic:.6f}  p = {pearson.pvalue:.6g}")
print(f"with Yates' correction, FOR COMPARISON ONLY: {yates.statistic:.6f}  "
      f"p = {yates.pvalue:.6g}")
print(f"Fisher's exact test: p = {fisher.pvalue:.6g}")

# THE ODDS RATIO IS COMPUTED FROM THE TABLE, NOT TAKEN FROM fisher_exact.
# scipy's fisher_exact returns the SAMPLE odds ratio ad/bc, while R's fisher.test returns the
# CONDITIONAL MAXIMUM LIKELIHOOD estimate -- a different estimator with a different value.
# scipy.stats.contingency.odds_ratio is the one that matches R's. Computing it from the table in
# both files removes the question.
log_or = float(np.log((a * e) / (b * c)))
log_or_se = float(np.sqrt(1 / a + 1 / b + 1 / c + 1 / e))
print(f"\nsample odds ratio {np.exp(log_or):.6f} "
      f"(95% CI {np.exp(log_or - 1.96 * log_or_se):.6f} to "
      f"{np.exp(log_or + 1.96 * log_or_se):.6f})")
print(f"fisher_exact's statistic is the SAME sample odds ratio: {fisher.statistic:.6f}")

print("\n--- HARNESS ---")
print(f"chisq_uncorrected={pearson.statistic:.10f}")
print(f"chisq_yates={yates.statistic:.10f}")
print(f"log_or={log_or:.10f}")
print(f"log_or_se={log_or_se:.10f}")
print(f"neg_log10_p_fisher={-np.log10(fisher.pvalue):.10f}")
print(f"risk_exposed={a / (a + b):.10f}")
print(f"risk_unexposed={c / (c + e):.10f}")
print(f"risk_difference={a / (a + b) - c / (c + e):.10f}")
print(f"n={len(d)}")
