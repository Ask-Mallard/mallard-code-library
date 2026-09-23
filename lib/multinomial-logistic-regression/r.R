# Multinomial logistic regression for an unordered categorical outcome
#
# One log relative-risk ratio per non-reference category: how treatment shifts the odds of rehab
# versus home, and of nursing versus home. The reference category decides what every coefficient
# means, so it is chosen deliberately and named in the report.

library(nnet)

d <- read.csv("fixture.csv", stringsAsFactors = FALSE)
stopifnot(!anyNA(d))

# PINNED: home is the reference. factor() on this text would put "home" first only because it sorts
# first alphabetically; an outcome like c("rehab", "home", "nursing") would silently make "home" the
# second level and change every coefficient. The levels are written out.
d$destination <- factor(d$destination, levels = c("home", "rehab", "nursing"))

# PINNED: maxit and reltol. multinom stops at 100 iterations and a loose relative tolerance by
# default, which here left the coefficients short of the maximum.
fit <- multinom(destination ~ treated, data = d, maxit = 1000, reltol = 1e-14, trace = FALSE)
stopifnot(fit$convergence == 0)

b <- coef(fit)
se <- summary(fit)$standard.errors
z <- qnorm(0.975)
print(b); print(se)

cat("\n--- HARNESS ---\n")
for (k in c("rehab", "nursing")) {
  cat(sprintf("%s_intercept=%.10f\n", k, b[k, "(Intercept)"]))
  cat(sprintf("%s_log_rrr=%.10f\n%s_log_rrr_se=%.10f\n", k, b[k, "treated"], k, se[k, "treated"]))
  cat(sprintf("%s_log_rrr_lcl=%.10f\n%s_log_rrr_ucl=%.10f\n", k, b[k, "treated"] - z * se[k, "treated"],
              k, b[k, "treated"] + z * se[k, "treated"]))
}
cat(sprintf("n=%d\n", nrow(d)))
