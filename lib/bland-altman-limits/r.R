# Method comparison: Bland-Altman bias and 95% limits of agreement, with their intervals.
#
# 1. CORRELATION IS NOT AGREEMENT. Two methods can correlate at 0.95 while one reads consistently higher
#    or the two differ by clinically important amounts. The questions are: how far apart are they on
#    average (bias), and within what range do individual differences fall (limits of agreement)?
#
# 2. THE LIMITS are bias +/- 1.96 SD of the differences. Whether they are narrow enough is a clinical
#    judgement made before the study; state the acceptable range in the protocol.
#
# 3. INTERVALS (Bland and Altman 1999): the bias with a t interval; each limit with SE = sqrt(3 s^2 / n)
#    and a t quantile on n - 1 df. Pinned: the approximate form, written out; the exact form (Carkeet 2015)
#    is slightly wider in small samples.
#
# 4. CHECK FOR PROPORTIONAL BIAS by regressing the differences on the means; a slope far from 0 means the
#    limits are not constant across the range, and a plot of differences against means should be shown.

d <- read.csv("fixture.csv")
dif <- d$method_b - d$method_a
avg <- (d$method_a + d$method_b) / 2
n <- length(dif)
bias <- mean(dif)
s <- sd(dif)
tq <- qt(0.975, n - 1)
zq <- qnorm(0.975)
loa <- bias + c(-1, 1) * zq * s
se_bias <- s / sqrt(n)
se_loa <- sqrt(3 * s^2 / n)
slope <- coef(lm(dif ~ avg))[["avg"]]

cat(sprintf("bias %.2f (95%% CI %.2f to %.2f); limits of agreement %.2f to %.2f; correlation %.3f\n",
            bias, bias - tq * se_bias, bias + tq * se_bias, loa[1], loa[2], cor(d$method_a, d$method_b)))

cat("\n--- HARNESS ---\n")
cat(sprintf("bias=%.10f\nbias_lcl=%.10f\nbias_ucl=%.10f\n", bias, bias - tq * se_bias, bias + tq * se_bias))
cat(sprintf("sd_difference=%.10f\n", s))
cat(sprintf("loa_lower=%.10f\nloa_lower_lcl=%.10f\nloa_lower_ucl=%.10f\n", loa[1], loa[1] - tq * se_loa, loa[1] + tq * se_loa))
cat(sprintf("loa_upper=%.10f\nloa_upper_lcl=%.10f\nloa_upper_ucl=%.10f\n", loa[2], loa[2] - tq * se_loa, loa[2] + tq * se_loa))
cat(sprintf("proportional_bias_slope=%.10f\ncorrelation=%.10f\nn=%d\n", slope, cor(d$method_a, d$method_b), n))
