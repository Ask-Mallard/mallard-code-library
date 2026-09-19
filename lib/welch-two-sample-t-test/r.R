# Two-sample comparison of means with unequal variances (Welch)
#
# Base R only. t.test is in stats, so this file installs nothing, which is worth noticing: the
# comparison clinicians run most often needs no package at all, and the thing that varies between
# languages is not capability but which test each one runs when you do not say.

d <- read.csv("fixture.csv")

# The arms are selected by the arm COLUMN, not by position. The fixture is shuffled so that an
# implementation reading "the first 1200 rows" fails visibly rather than by luck.
treated <- d$y[d$arm == 1]
control <- d$y[d$arm == 0]

# PINNED DEFAULT: var.equal = FALSE, the Welch test.
#
# This IS R's default and it is named anyway, because it is NOT scipy's default:
#   R      t.test()               Welch
#   Python scipy.stats.ttest_ind  POOLED -- equal_var=True is the default
#
# On this fixture the pooled standard error is 32% smaller than Welch's, because the smaller arm is
# the more variable one. That direction makes the pooled test anti-conservative: it manufactures
# significance rather than losing it.
tt <- t.test(treated, control, var.equal = FALSE)

# The pooled test, computed only so the gap can be reported. It is not the recommended analysis
# here and no interval in this file is built from it.
pooled <- t.test(treated, control, var.equal = TRUE)

cat("Welch (unequal variances)\n")
print(tt)
cat("\npooled, FOR COMPARISON ONLY -- it assumes the arms have the same spread\n")
print(pooled)

cat("\n--- HARNESS ---\n")
cat(sprintf("mean_difference=%.10f\n", mean(treated) - mean(control)))
cat(sprintf("mean_treated=%.10f\n", mean(treated)))
cat(sprintf("mean_control=%.10f\n", mean(control)))
cat(sprintf("difference_se=%.10f\n", tt$stderr))
cat(sprintf("difference_se_pooled=%.10f\n", pooled$stderr))
cat(sprintf("welch_df=%.10f\n", unname(tt$parameter)))
cat(sprintf("pooled_df=%.10f\n", unname(pooled$parameter)))
cat(sprintf("t_welch=%.10f\n", unname(tt$statistic)))
cat(sprintf("n_treated=%d\n", length(treated)))
cat(sprintf("n_control=%d\n", length(control)))
