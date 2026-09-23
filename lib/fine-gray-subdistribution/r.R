# Fine-Gray regression: the effect of treatment on the CUMULATIVE INCIDENCE of one cause
#
# Under competing risks, the question "does treatment change how many people relapse?" is about the
# cumulative incidence, and the Fine-Gray model answers it with a subdistribution hazard ratio. Its
# risk set keeps people who already died of the competing cause, which is why the ratio is not a rate
# among those at risk; interpret it through the cumulative incidence it implies. For the effect on
# the RATE of relapse among those still alive, see cause-specific-cox. Reporting both is common.

suppressPackageStartupMessages(library(cmprsk))

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d), all(d$status %in% c(0, 1, 2)))

# PINNED: failcode = 1 (the cause of interest) and cencode = 0, named although they are crr's defaults,
# because a data set coding death as 1 would silently model death. cov1 must be a numeric MATRIX:
# crr does not take a formula or expand factors.
fit <- crr(ftime = d$time, fstatus = d$status, cov1 = cbind(treated = d$treated),
           failcode = 1, cencode = 0)
b <- fit$coef[["treated"]]; se <- sqrt(fit$var[1, 1])
z <- qnorm(0.975)
cat(sprintf("subdistribution HR for relapse %.3f, 95%% %.3f to %.3f\n", exp(b), exp(b - z * se), exp(b + z * se)))

cat("\n--- HARNESS ---\n")
cat(sprintf("log_subdistribution_hr=%.10f\nlog_subdistribution_hr_se=%.10f\n", b, se))
cat(sprintf("log_subdistribution_hr_lcl=%.10f\nlog_subdistribution_hr_ucl=%.10f\n", b - z * se, b + z * se))
cat(sprintf("cause1=%d\ncause2=%d\nn=%d\n", sum(d$status == 1), sum(d$status == 2), nrow(d)))
