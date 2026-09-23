# Two-period crossover (AB/BA): the period-adjusted treatment effect
#
# Each patient's period-1 minus period-2 difference contains the treatment effect with opposite signs
# in the two sequences and the period effect with the same sign. Comparing those differences between
# sequences cancels the period effect: treatment effect (B - A) = (mean diff in BA - mean diff in AB) / 2.
# It assumes no carryover, which the washout must make plausible; the data cannot test it well.

d <- read.csv("fixture.csv", stringsAsFactors = FALSE)
stopifnot(!anyNA(d))

w <- reshape(d[, c("id", "sequence", "period", "y")], idvar = c("id", "sequence"), timevar = "period",
             direction = "wide")
w$diff <- w$y.1 - w$y.2

# PINNED: the two-sample t-test of the period differences BY SEQUENCE (BA first), halved. A paired
# t-test of B against A ignores the period, and with unequal sequences the period effect biases it.
fit <- t.test(diff ~ factor(sequence, levels = c("BA", "AB")), data = w, var.equal = TRUE, conf.level = 0.95)
effect <- unname(fit$estimate[1] - fit$estimate[2]) / 2
ci <- fit$conf.int / 2

cat(sprintf("treatment effect (B - A) %.3f, 95%% %.3f to %.3f, %d df\n", effect, ci[1], ci[2], as.integer(fit$parameter)))

cat("\n--- HARNESS ---\n")
cat(sprintf("treatment_effect=%.10f\ntreatment_effect_lcl=%.10f\ntreatment_effect_ucl=%.10f\n", effect, ci[1], ci[2]))
cat(sprintf("df=%d\npatients=%d\n", as.integer(fit$parameter), nrow(w)))
