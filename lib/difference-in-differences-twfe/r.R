# Difference-in-differences: two-way fixed effects for a policy adopted at one time.
#
# 1. THE COMPARISON IS OF CHANGES, NOT LEVELS. The adopters' before-after change contains the secular
#    trend; the after-period difference between groups contains their baseline difference. The DiD
#    subtracts the non-adopters' change from the adopters' change, which removes both, provided the two
#    groups would have followed PARALLEL TRENDS without the policy. That is an assumption; plot the
#    pre-period means of both groups (or fit an event study) before believing it.
#
# 2. TWO-WAY FIXED EFFECTS, HOSPITAL AND QUARTER. With one adoption date and a balanced panel the
#    coefficient on the policy indicator equals the 2x2 DiD of group means exactly (the control file
#    checks it). With STAGGERED adoption and effects that change over time it does not, and can even
#    take the wrong sign (Goodman-Bacon 2021); use a heterogeneity-robust estimator there
#    (Callaway-Sant'Anna, did package). This entry is the single-date case only.
#
# 3. CLUSTER BY THE UNIT THAT ADOPTS. Quarters within a hospital are correlated; the SE clusters on
#    hospital (sandwich::vcovCL, HC1). 40 clusters here; with fewer, see the CR2 entry.

library(sandwich)

d <- read.csv("fixture.csv")
fit <- lm(los ~ policy + factor(hospital) + factor(quarter), data = d)
did <- coef(fit)[["policy"]]
se <- sqrt(vcovCL(fit, cluster = ~hospital, type = "HC1")["policy", "policy"])

# The two comparisons DiD exists to replace.
g <- function(a, p) mean(d$los[d$adopter == a & d$post == p])
pre_post_adopters <- g(1, 1) - g(1, 0)
post_between_groups <- g(1, 1) - g(0, 1)

z <- qnorm(0.975)
cat(sprintf("DiD %.3f (95%% CI %.3f to %.3f)\n", did, did - z * se, did + z * se))

cat("\n--- HARNESS ---\n")
cat(sprintf("did=%.10f\ndid_se=%.10f\ndid_lcl=%.10f\ndid_ucl=%.10f\n", did, se, did - z * se, did + z * se))
cat(sprintf("pre_post_adopters=%.10f\npost_between_groups=%.10f\n", pre_post_adopters, post_between_groups))
cat(sprintf("hospitals=%d\nn=%d\n", length(unique(d$hospital)), nrow(d)))
