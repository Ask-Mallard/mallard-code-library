# A single proportion (a prevalence) with Wilson and Clopper-Pearson 95% intervals
#
# This is the one case where these intervals belong. They are for a SINGLE proportion or rate from
# independent observations. An odds ratio, a risk ratio or an adjusted model estimate takes a
# profile-likelihood or Wald interval instead, and clustered data takes a design-based interval
# (see clustered-prevalence-cluster-sample).
#
# Wilson is the usual report. Clopper-Pearson ("exact") guarantees at least nominal coverage and is
# wider; report it when a conservative interval is wanted, for example with very few events.

library(binom)

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d), all(d$outcome %in% c(0, 1)))

x <- sum(d$outcome)
n <- nrow(d)

# PINNED DEFAULT: methods = "wilson". binom.confint's default is methods = "all", which returns
# eleven intervals in a data frame, and a caller who indexes the first row reports whichever method
# happens to sort first.
wilson <- binom.confint(x, n, conf.level = 0.95, methods = "wilson")

# Clopper-Pearson from base R's binom.test, an implementation independent of the binom package.
cp <- binom.test(x, n, conf.level = 0.95)$conf.int

cat(sprintf("events %d of %d, proportion %.4f\n", x, n, x / n))
cat(sprintf("Wilson 95%%          %.4f to %.4f\n", wilson$lower, wilson$upper))
cat(sprintf("Clopper-Pearson 95%% %.4f to %.4f\n", cp[1], cp[2]))

cat("\n--- HARNESS ---\n")
cat(sprintf("prevalence=%.10f\n", x / n))
cat(sprintf("wilson_lcl=%.10f\n", wilson$lower))
cat(sprintf("wilson_ucl=%.10f\n", wilson$upper))
cat(sprintf("cp_lcl=%.10f\n", cp[1]))
cat(sprintf("cp_ucl=%.10f\n", cp[2]))
cat(sprintf("events=%d\n", x))
cat(sprintf("n=%d\n", n))
