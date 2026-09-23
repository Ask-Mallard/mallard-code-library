# Restricted mean survival time (RMST) and the difference between arms, up to a stated horizon
#
# The area under each Kaplan-Meier curve up to tau: the mean event-free time over those tau years. The
# difference ("treated patients lived 0.42 years longer on average over 4 years") needs no
# proportional-hazards assumption and is in units a patient understands.

library(survival)

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d), all(d$arm %in% c(0, 1)))

# PINNED: rmean = 4, the horizon stated before looking at the data. print.survfit's default reports
# no restricted mean, and "rmean = 'common'" picks the largest observed time, which moves with the
# data. tau must lie within follow-up in BOTH arms.
tau <- 4
tab <- summary(survfit(Surv(time, event) ~ arm, data = d), rmean = tau)$table
r1 <- tab["arm=1", "rmean"]; s1 <- tab["arm=1", "se(rmean)"]
r0 <- tab["arm=0", "rmean"]; s0 <- tab["arm=0", "se(rmean)"]
diff <- r1 - r0; se <- sqrt(s1^2 + s0^2); z <- qnorm(0.975)

cat(sprintf("RMST to %g years: treated %.3f, control %.3f; difference %.3f, 95%% %.3f to %.3f\n",
            tau, r1, r0, diff, diff - z * se, diff + z * se))

cat("\n--- HARNESS ---\n")
cat(sprintf("rmst_treated=%.10f\nrmst_treated_se=%.10f\n", r1, s1))
cat(sprintf("rmst_control=%.10f\nrmst_control_se=%.10f\n", r0, s0))
cat(sprintf("rmst_difference=%.10f\nrmst_difference_se=%.10f\n", diff, se))
cat(sprintf("rmst_difference_lcl=%.10f\nrmst_difference_ucl=%.10f\n", diff - z * se, diff + z * se))
cat(sprintf("n=%d\n", nrow(d)))
