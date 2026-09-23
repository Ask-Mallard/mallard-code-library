# An incidence rate (events per person-year) with an exact Poisson 95% interval
#
# The denominator is PERSON-TIME, not people. A cohort where people are followed for different
# lengths of time has no single "risk"; the rate is what it estimates. The exact interval assumes
# the events are Poisson: independent, at one rate for everyone. When some people have repeated
# events, the data are overdispersed and this interval is too narrow (see README.md).

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d), all(d$followup_years > 0), all(d$events >= 0))

events <- sum(d$events)
person_years <- sum(d$followup_years)

# PINNED DEFAULT: T = person_years. poisson.test's default is T = 1, so omitting it returns a
# perfectly valid interval for the COUNT of events, not for the rate, with no warning.
fit <- poisson.test(events, T = person_years, conf.level = 0.95)

cat(sprintf("%d events over %.1f person-years\n", events, person_years))
cat(sprintf("rate %.4f per person-year (%.1f per 1000), exact 95%% %.4f to %.4f\n",
            fit$estimate, 1000 * fit$estimate, fit$conf.int[1], fit$conf.int[2]))

cat("\n--- HARNESS ---\n")
cat(sprintf("rate=%.10f\n", fit$estimate))
cat(sprintf("rate_lcl=%.10f\n", fit$conf.int[1]))
cat(sprintf("rate_ucl=%.10f\n", fit$conf.int[2]))
cat(sprintf("events=%d\n", events))
cat(sprintf("person_years=%.10f\n", person_years))
cat(sprintf("n=%d\n", nrow(d)))
