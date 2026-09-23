# Proportional-odds (cumulative logit) regression for an ordered outcome
#
# The Python equivalent of r.R, using statsmodels' OrderedModel, independent of R's polr. Both use
# P(Y <= k) = F(threshold_k - x'beta), so a positive coefficient means higher categories.

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.miscmodels.ordinal_model import OrderedModel

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any()

# PINNED: the category order written out; pandas would otherwise order the text alphabetically.
levels = ["none", "mild", "moderate", "severe"]
y = pd.Series(pd.Categorical(d.severity, categories=levels, ordered=True))
assert not y.isna().any()

model = OrderedModel(y, d[["treated"]].astype(float), distr="logit")
fit = model.fit(method="newton", disp=False, maxiter=200, tol=1e-12)
assert fit.mle_retvals["converged"], "the ordinal fit did not converge"

z = stats.norm.ppf(0.975)
b, se = fit.params["treated"], fit.bse["treated"]
# The thresholds are stored as the first cut plus log-increments; transform them back.
cuts = model.transform_threshold_params(fit.params.to_numpy())[1:-1]
print(f"odds ratio for a higher category {np.exp(b):.3f}, Wald 95% {np.exp(b - z * se):.3f} to {np.exp(b + z * se):.3f}")
print(cuts)

print("\n--- HARNESS ---")
print(f"log_odds_ratio={b:.10f}\nlog_odds_ratio_se={se:.10f}")
print(f"log_odds_ratio_lcl={b - z * se:.10f}\nlog_odds_ratio_ucl={b + z * se:.10f}")
for k, c in enumerate(cuts, 1):
    print(f"threshold_{k}={c:.10f}")
print(f"n={len(d)}")
