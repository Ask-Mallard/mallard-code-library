# Marginal logistic regression by GEE, on binary outcomes clustered within clinic.
#
# THREE THINGS ARE PINNED HERE AND EVERY ONE OF THEM IS A DEFAULT SOMEWHERE ELSE.
#
# 1. corstr = "exchangeable". geeglm DEFAULTS TO "independence", and so does statsmodels. So the
#    same analysis fits a different working correlation depending on whether the default is written
#    out, and no output says which.
#
# 2. The data are SORTED BY CLUSTER. geeglm requires observations from one cluster to be
#    contiguous; it does not check, and unsorted input is silently treated as many small clusters.
#    That is the failure mode this file guards against by sorting unconditionally, even though the
#    committed fixture already arrives sorted.
#
# 3. BOTH standard errors are reported, under BOTH working correlations. geeglm's summary prints
#    the sandwich; the model-based error is the alternative. How far apart those two are DEPENDS ON
#    THE WORKING CORRELATION, which is why this file reports four numbers and not two: a model-based
#    variance computed under an exchangeable structure already carries the within-clinic correlation
#    and lands near the sandwich, while one computed under independence does not. The dangerous
#    combination is independence with a model-based error -- a pairing nobody would write on purpose,
#    but one a careless edit produces by taking the model-based error under the default independence
#    structure.

library(geepack)

d <- read.csv("fixture.csv")
d <- d[order(d$clinic), ]

fit <- geeglm(outcome ~ exposed + x, id = clinic, data = d,
              family = binomial(link = "logit"), corstr = "exchangeable")

# vbeta is the sandwich, vbeta.naiv the model-based variance. Taken from the fitted object rather
# than from summary(), because summary() prints one of them and the entry needs both.
se_robust <- sqrt(diag(fit$geese$vbeta))
se_naive <- sqrt(diag(fit$geese$vbeta.naiv))
nm <- names(coef(fit))

# The independence working correlation, which is what geeglm, statsmodels and PROC GENMOD all
# give you by default. Both of ITS standard errors are reported too, because the size of the
# sandwich-to-model-based gap depends on which working correlation produced it: a model-based
# variance that already carries the correlation is close to the sandwich, and one that assumes
# independence when the data are correlated is not.
ind <- geeglm(outcome ~ exposed + x, id = clinic, data = d,
              family = binomial(link = "logit"), corstr = "independence")
ind_robust <- sqrt(diag(ind$geese$vbeta))
ind_naive <- sqrt(diag(ind$geese$vbeta.naiv))
ind_nm <- names(coef(ind))

cat("--- HARNESS ---\n")
cat(sprintf("exposure_log_or=%.10f\n", coef(fit)[["exposed"]]))
cat(sprintf("exposure_se_robust=%.10f\n", se_robust[which(nm == "exposed")]))
cat(sprintf("exposure_se_naive=%.10f\n", se_naive[which(nm == "exposed")]))
cat(sprintf("covariate_beta=%.10f\n", coef(fit)[["x"]]))
cat(sprintf("exposure_log_or_independence=%.10f\n", coef(ind)[["exposed"]]))
cat(sprintf("exposure_se_robust_independence=%.10f\n", ind_robust[which(ind_nm == "exposed")]))
cat(sprintf("exposure_se_naive_independence=%.10f\n", ind_naive[which(ind_nm == "exposed")]))
cat(sprintf("alpha_exchangeable=%.10f\n", fit$geese$alpha))
