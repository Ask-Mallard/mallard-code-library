# Median (quantile) regression: the adjusted difference in MEDIAN length of stay
#
# For a skewed outcome, the median effect is often the question a clinician means ("a typical patient
# stays 1.5 days longer"), and it is a different number from the mean effect when the spread differs
# between groups. Median regression estimates it with covariate adjustment, which a rank test cannot.

suppressPackageStartupMessages(library(quantreg))

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d))
d$age60 <- d$age - 60

# PINNED: tau = 0.5 and the Barrodale-Roberts simplex ("br", exact). PINNED: se = "nid", the
# Hendricks-Koenker local sandwich with the Hall-Sheather bandwidth. summary.rq's default for n < 1001
# is "rank" (inverting a rank test), which gives an interval but no standard error; "iid" assumes one
# error distribution for everybody, which is false here, where the treated vary more.
fit <- rq(los_days ~ treated + age60, tau = 0.5, data = d, method = "br")
s <- summary(fit, se = "nid", hs = TRUE)$coefficients
tq <- qt(0.975, nrow(d) - length(coef(fit)))
print(s)

cat("\n--- HARNESS ---\n")
cat(sprintf("median_effect=%.10f\nmedian_effect_se=%.10f\n", s["treated", "Value"], s["treated", "Std. Error"]))
cat(sprintf("median_effect_lcl=%.10f\nmedian_effect_ucl=%.10f\n",
            s["treated", "Value"] - tq * s["treated", "Std. Error"], s["treated", "Value"] + tq * s["treated", "Std. Error"]))
cat(sprintf("age_slope=%.10f\n", s["age60", "Value"]))
cat(sprintf("n=%d\n", nrow(d)))
