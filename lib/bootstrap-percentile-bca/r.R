# Bootstrap confidence intervals for a ratio of means from skewed data: percentile and BCa.
#
# 1. RESAMPLE AS THE DATA WERE SAMPLED. Two independent hospitals: resample patients within each hospital
#    (a stratified bootstrap), recompute the ratio of mean length of stay, repeat B = 2,000 times.
#
# 2. THE INTERVAL.
#    - Percentile: the 2.5th and 97.5th percentiles of the bootstrap distribution. Simple; biased when the
#      statistic is biased or its spread changes with its value.
#    - BCa (bias-corrected and accelerated; Efron 1987): shifts the percentiles by a bias correction z0
#      (the share of bootstrap values below the estimate) and an acceleration a (from the jackknife). It is
#      second-order accurate, and the one to report by default. Its jackknife acceleration needs a smooth
#      statistic: for a median it is exactly 0 (see fixture.py). boot::boot.ci(type = "bca") computes it.
#
# 3. THE RESAMPLES HERE come from a shared "minimal standard" generator (Park and Miller, multiplier 48271,
#    modulus 2^31 - 1), implemented identically in the Python file, so both languages use the same
#    resamples and agree exactly. In practice use R's own generator (boot::boot with a set seed).
#    Quantiles are R's type 7 (numpy's default).

d <- read.csv("fixture.csv")
a <- d$los[d$hospital == "A"]
b <- d$los[d$hospital == "B"]
estimate <- mean(a) / mean(b)
B <- 2000L

state <- 20261201
draw <- function(n) {  # n uniform indices in 1..n from the shared generator
  idx <- integer(n)
  for (k in seq_len(n)) {
    state <<- (48271 * state) %% 2147483647
    idx[k] <- floor(state / 2147483647 * n) + 1L
  }
  idx
}
boot_stat <- numeric(B)
for (r in seq_len(B)) boot_stat[r] <- mean(a[draw(length(a))]) / mean(b[draw(length(b))])

percentile <- quantile(boot_stat, c(0.025, 0.975), type = 7, names = FALSE)

# BCa: bias correction and jackknife acceleration over all patients.
z0 <- qnorm(mean(boot_stat < estimate))
jack <- c(sapply(seq_along(a), function(i) mean(a[-i]) / mean(b)),
          sapply(seq_along(b), function(i) mean(a) / mean(b[-i])))
jm <- mean(jack)
acc <- sum((jm - jack)^3) / (6 * sum((jm - jack)^2)^1.5)
alpha <- c(0.025, 0.975)
adj <- pnorm(z0 + (z0 + qnorm(alpha)) / (1 - acc * (z0 + qnorm(alpha))))
bca <- quantile(boot_stat, adj, type = 7, names = FALSE)

cat(sprintf("ratio of means %.3f; percentile 95%% CI %.3f to %.3f; BCa %.3f to %.3f; bootstrap SE %.3f\n",
            estimate, percentile[1], percentile[2], bca[1], bca[2], sd(boot_stat)))

cat("\n--- HARNESS ---\n")
cat(sprintf("mean_ratio=%.10f\nbootstrap_se=%.10f\n", estimate, sd(boot_stat)))
cat(sprintf("percentile_lcl=%.10f\npercentile_ucl=%.10f\n", percentile[1], percentile[2]))
cat(sprintf("bca_lcl=%.10f\nbca_ucl=%.10f\nbias_correction_z0=%.10f\nacceleration=%.10f\n", bca[1], bca[2], z0, acc))
cat(sprintf("resamples=%d\n", B))
