# Network meta-analysis (frequentist, random effects) of three treatments from two-arm trials.
#
# 1. A NETWORK COMBINES DIRECT AND INDIRECT EVIDENCE. Trials of A vs B and A vs C also tell us about C vs
#    B, through A. The combined estimate is valid only if the network is CONSISTENT (direct and indirect
#    evidence agree) and the trials are similar enough in effect modifiers (transitivity), which is a
#    clinical judgement. Report the test of inconsistency beside the estimates.
#
# 2. THE MODEL is netmeta's (Rucker 2012): weighted least squares on the trial contrasts, a common
#    between-trial variance tau^2 by the generalised DerSimonian-Laird estimator, random effects. For
#    two-arm trials it is written out here in a few lines, and reproduces netmeta::netmeta(sm = "OR",
#    method.tau = "DL", reference.group = "A") to 1e-10 on this fixture (checked with netmeta 3.7.0).
#    netmeta itself is not called: it depends on igraph, whose compiled library needs a system library
#    the CI runner does not have. With multi-arm trials, use netmeta, which handles their correlation.
#
# 3. RANKINGS (P-scores, SUCRA) summarize uncertainty poorly on their own; report effects with intervals
#    first.

d <- read.csv("fixture.csv")
params <- c("B", "C")  # basic parameters against the reference A
X <- matrix(0, nrow(d), 2, dimnames = list(NULL, params))
for (i in seq_len(nrow(d))) {
  if (d$treat1[i] %in% params) X[i, d$treat1[i]] <- X[i, d$treat1[i]] + 1
  if (d$treat2[i] %in% params) X[i, d$treat2[i]] <- X[i, d$treat2[i]] - 1
}
y <- d$TE
s2 <- d$seTE^2

wls <- function(w) {
  info <- crossprod(X * w, X)
  list(beta = drop(solve(info, crossprod(X * w, y))), cov = solve(info))
}
w <- 1 / s2
common <- wls(w)
q <- sum(w * (y - X %*% common$beta)^2)
df <- length(y) - ncol(X)
H <- (X * w) %*% solve(crossprod(X * w, X)) %*% t(X * w)
tau2 <- max(0, (q - df) / (sum(w) - sum(diag(H))))  # PINNED: generalised DerSimonian-Laird
re <- wls(1 / (s2 + tau2))

# Within-design heterogeneity; the rest of Q is between designs (inconsistency).
design <- paste(d$treat1, d$treat2, sep = ":")
q_within <- sum(sapply(split(seq_along(y), design), function(i) {
  mu <- sum(w[i] * y[i]) / sum(w[i]); sum(w[i] * (y[i] - mu)^2)
}))

contrasts <- list(B_vs_A = c(1, 0), C_vs_A = c(0, 1), C_vs_B = c(-1, 1))
est <- sapply(contrasts, function(c) sum(c * re$beta))
se <- sapply(contrasts, function(c) sqrt(drop(t(c) %*% re$cov %*% c)))
z <- qnorm(0.975)

for (k in c("B_vs_A", "C_vs_A")) cat(sprintf("%s: OR %.3f (95%% CI %.3f to %.3f)\n", sub("_vs_", " vs ", k),
                                             exp(est[[k]]), exp(est[[k]] - z * se[[k]]), exp(est[[k]] + z * se[[k]])))
cat(sprintf("C vs B: OR %.3f; tau^2 %.4f; Q %.3f (between designs %.3f)\n", exp(est[["C_vs_B"]]), tau2, q, q - q_within))

cat("\n--- HARNESS ---\n")
for (k in names(contrasts)) cat(sprintf("log_or_%s=%.10f\nlog_or_%s_se=%.10f\n", k, est[[k]], k, se[[k]]))
cat(sprintf("common_log_or_C_vs_B=%.10f\n", sum(c(-1, 1) * common$beta)))
cat(sprintf("tau2=%.10f\nq=%.10f\nq_inconsistency=%.10f\n", tau2, q, q - q_within))
cat(sprintf("studies=%d\n", nrow(d)))
