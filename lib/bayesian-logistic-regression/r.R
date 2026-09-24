# Bayesian logistic regression with stated priors, sampled by Markov chain Monte Carlo.
#
# 1. STATE THE PRIORS, and why. Here: normal(0, 2.5) on each coefficient (weakly informative on the
#    log-odds scale: it rules out odds ratios beyond about e^5 without favouring a direction) and
#    normal(0, 5) on the intercept. Report a sensitivity analysis under a sceptical prior where the
#    conclusion could depend on it.
#
# 2. REPORT THE POSTERIOR, not a p-value: the posterior mean and SD of the treatment log odds ratio, a
#    95% credible interval, and the posterior probability that the effect is beneficial or harmful.
#
# 3. CHECK THE SAMPLER. Several chains from dispersed starts; R-hat (Gelman-Rubin) near 1 and an effective
#    sample size in the thousands before reading the numbers.
#
# 4. THE SAMPLER here is random-walk Metropolis written out (4 chains of 5,000 after 1,000 warm-up; the
#    proposal covariance is 2.4^2 / p times the maximum-likelihood covariance), so the entry needs no
#    compiled software and every step is visible. In practice use Stan (rstanarm::stan_glm or brms),
#    whose sampler is far more efficient; the model and priors are the same.
#
# R ONLY, TRUTH CHECK ONLY: Markov chains cannot be matched draw for draw across
# languages. The control file checks the posterior against a Laplace approximation and recovery over seeds.

set.seed(20261206)
d <- read.csv("fixture.csv")
X <- cbind(1, d$treated, d$age)
y <- d$event
prior_sd <- c(5, 2.5, 2.5)  # PINNED: intercept, treatment, age

log_post <- function(b) {
  eta <- drop(X %*% b)
  sum(y * eta - log1p(exp(eta))) + sum(dnorm(b, 0, prior_sd, log = TRUE))
}

mle <- glm(event ~ treated + age, data = d, family = binomial)
prop_chol <- chol(2.4^2 / ncol(X) * vcov(mle))
chains <- 4; warmup <- 1000; keep <- 5000
draws <- array(NA_real_, c(keep, chains, ncol(X)))
accepted <- 0
for (ch in seq_len(chains)) {
  b <- coef(mle) + rnorm(ncol(X), 0, 2 * sqrt(diag(vcov(mle))))  # dispersed start
  lp <- log_post(b)
  for (it in seq_len(warmup + keep)) {
    cand <- b + drop(rnorm(ncol(X)) %*% prop_chol)
    lc <- log_post(cand)
    if (log(runif(1)) < lc - lp) { b <- cand; lp <- lc; if (it > warmup) accepted <- accepted + 1 }
    if (it > warmup) draws[it - warmup, ch, ] <- b
  }
}

theta <- draws[, , 2]
# The classic Gelman-Rubin R-hat over the 4 chains, and an effective sample size from the chains'
# average autocorrelation, summed until it falls below 0.05.
W <- mean(apply(theta, 2, var)); Bv <- keep * var(colMeans(theta))
rhat <- sqrt(((keep - 1) / keep * W + Bv / keep) / W)
rho <- rowMeans(sapply(seq_len(chains), function(ch) acf(theta[, ch], lag.max = 200, plot = FALSE)$acf[-1]))
ess <- length(theta) / (1 + 2 * sum(rho[seq_len(which(rho < 0.05)[1] - 1)]))
post <- as.vector(theta)

cat(sprintf("treatment log OR: posterior mean %.3f (SD %.3f), 95%% credible interval %.3f to %.3f; P(log OR > 0) = %.3f; R-hat %.3f, ESS about %.0f\n",
            mean(post), sd(post), quantile(post, 0.025), quantile(post, 0.975), mean(post > 0), rhat, ess))

cat("\n--- HARNESS ---\n")
cat(sprintf("treatment_log_or=%.10f\nposterior_sd=%.10f\n", mean(post), sd(post)))
cat(sprintf("credible_lcl=%.10f\ncredible_ucl=%.10f\nprob_positive=%.10f\n", quantile(post, 0.025), quantile(post, 0.975), mean(post > 0)))
cat(sprintf("rhat=%.10f\ness=%.10f\nacceptance=%.10f\n", rhat, ess, accepted / (chains * keep)))
cat(sprintf("mle_treatment_log_or=%.10f\nn=%d\n", coef(mle)[["treated"]], nrow(d)))
