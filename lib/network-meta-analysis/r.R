# Network meta-analysis (frequentist, random effects) of three treatments from two-arm trials.
#
# 1. A NETWORK COMBINES DIRECT AND INDIRECT EVIDENCE. Trials of A vs B and A vs C also tell us about C vs
#    B, through A. The combined estimate is valid only if the network is CONSISTENT (direct and indirect
#    evidence agree) and the trials are similar enough in effect modifiers (transitivity), which is a
#    clinical judgement. Report the test of inconsistency beside the estimates.
#
# 2. THE MODEL (netmeta, Rucker 2012): weighted least squares on the trial contrasts, a common
#    between-trial variance tau^2 by the generalised DerSimonian-Laird estimator, random effects.
#    Pinned: summary measure log OR, reference treatment A (placebo), random effects reported, the
#    common-effect Q and its between-design (inconsistency) component beside them.
#
# 3. RANKINGS (P-scores, SUCRA) summarize uncertainty poorly on their own; report effects with intervals
#    first.

suppressPackageStartupMessages(library(netmeta))

d <- read.csv("fixture.csv")
nm <- netmeta(TE, seTE, treat1, treat2, studlab = study, data = d, sm = "OR", reference.group = "A",
              common = TRUE, random = TRUE, method.tau = "DL")
re <- nm$TE.random
se <- nm$seTE.random
z <- qnorm(0.975)

for (tr in c("B", "C")) cat(sprintf("%s vs A: OR %.3f (95%% CI %.3f to %.3f)\n", tr, exp(re[tr, "A"]),
                                    exp(re[tr, "A"] - z * se[tr, "A"]), exp(re[tr, "A"] + z * se[tr, "A"])))
cat(sprintf("C vs B: OR %.3f; tau^2 %.4f; Q %.3f (between designs %.3f)\n", exp(re["C", "B"]), nm$tau^2, nm$Q, nm$Q.inconsistency))

cat("\n--- HARNESS ---\n")
cat(sprintf("log_or_B_vs_A=%.10f\nlog_or_B_vs_A_se=%.10f\n", re["B", "A"], se["B", "A"]))
cat(sprintf("log_or_C_vs_A=%.10f\nlog_or_C_vs_A_se=%.10f\n", re["C", "A"], se["C", "A"]))
cat(sprintf("log_or_C_vs_B=%.10f\nlog_or_C_vs_B_se=%.10f\n", re["C", "B"], se["C", "B"]))
cat(sprintf("common_log_or_C_vs_B=%.10f\n", nm$TE.common["C", "B"]))
cat(sprintf("tau2=%.10f\nq=%.10f\nq_inconsistency=%.10f\n", nm$tau^2, nm$Q, nm$Q.inconsistency))
cat(sprintf("studies=%d\n", nrow(d)))
