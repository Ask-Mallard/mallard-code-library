# Proportional-odds (cumulative logit) regression for an ordered outcome
#
# One odds ratio summarizes the shift toward HIGHER categories at every cut point together (none vs
# mild-or-worse, up to severe vs the rest). It assumes that odds ratio is the same at every cut: the
# proportional-odds assumption, which should be checked (e.g. by fitting the binary splits
# separately) and stated.

library(MASS)

d <- read.csv("fixture.csv", stringsAsFactors = FALSE)
stopifnot(!anyNA(d))

# PINNED: the level ORDER is written out. factor() on this text sorts alphabetically
# (mild, moderate, none, severe) and polr would fit that ordering without complaint.
d$severity <- factor(d$severity, levels = c("none", "mild", "moderate", "severe"), ordered = TRUE)
stopifnot(!anyNA(d$severity))

# polr's parameterization: logit P(Y <= k) = zeta_k - eta, so a POSITIVE coefficient means higher
# categories. Hess = TRUE so standard errors come from the fitted Hessian.
fit <- polr(severity ~ treated, data = d, method = "logistic", Hess = TRUE,
            control = list(reltol = 1e-14, maxit = 1000))

z <- qnorm(0.975)
b <- coef(fit)[["treated"]]
se <- sqrt(vcov(fit)["treated", "treated"])
cat(sprintf("odds ratio for a higher category %.3f, Wald 95%% %.3f to %.3f\n", exp(b), exp(b - z * se), exp(b + z * se)))
print(fit$zeta)

cat("\n--- HARNESS ---\n")
cat(sprintf("log_odds_ratio=%.10f\nlog_odds_ratio_se=%.10f\n", b, se))
cat(sprintf("log_odds_ratio_lcl=%.10f\nlog_odds_ratio_ucl=%.10f\n", b - z * se, b + z * se))
cat(sprintf("threshold_1=%.10f\nthreshold_2=%.10f\nthreshold_3=%.10f\n", fit$zeta[1], fit$zeta[2], fit$zeta[3]))
cat(sprintf("n=%d\n", nrow(d)))
