# Propensity-score matching: the average treatment effect in the TREATED.
#
# 1. MATCHING ESTIMATES THE ATT, NOT THE ATE. Each treated patient is kept and given a similar control,
#    so the answer describes the treated. When the effect varies with a confounder the two differ; here
#    the ATE is 2 and the ATT about 2.7. A plan that asks for "the effect of treatment" in everyone
#    wants weighting or g-computation, not this.
#
# 2. EVERY MATCHING OPTION IS NAMED. MatchIt's defaults are 1:1 nearest neighbour without replacement
#    on the propensity score, with no caliper. This file pins the distance as the LOGIT of the score
#    (link = "linear.logit"), a caliper of 0.2 SD of that logit (Austin 2011), the order in which treated
#    patients are matched (largest score first), and ATT. A different order or caliper gives a
#    different matched set and a different answer.
#
# 3. BALANCE IS THE DIAGNOSTIC, NOT A p-VALUE. The standardized mean difference of each covariate is
#    reported before and after matching, divided by the treated group's SD (the ATT convention).
#    Under 0.1 is the usual target.
#
# 4. THE SE COMES FROM THE PAIRS. A matched pair shares a propensity score; sandwich::vcovCL clustered
#    on the pair (subclass) is MatchIt's recommended variance. It treats the score as known, which is
#    conservative for the ATT (Abadie and Imbens 2016).

library(MatchIt)
library(sandwich)

d <- read.csv("fixture.csv")

m <- matchit(treated ~ severity + comorbid, data = d, method = "nearest", distance = "glm",
             link = "linear.logit", caliper = 0.2, std.caliper = TRUE, m.order = "largest",
             replace = FALSE, ratio = 1, estimand = "ATT")
md <- match.data(m)

fit <- lm(y ~ treated, data = md, weights = weights)
att <- coef(fit)[["treated"]]
se <- sqrt(vcovCL(fit, cluster = ~subclass)["treated", "treated"])

# Balance: SMD of severity with the treated group's SD from the full sample, before and after.
sd_t <- sd(d$severity[d$treated == 1])
smd <- function(x) (mean(x$severity[x$treated == 1]) - mean(x$severity[x$treated == 0])) / sd_t
naive <- mean(d$y[d$treated == 1]) - mean(d$y[d$treated == 0])

z <- qnorm(0.975)
cat(sprintf("ATT %.3f (95%% CI %.3f to %.3f); %d of %d treated matched\n", att, att - z * se, att + z * se,
            sum(md$treated == 1), sum(d$treated == 1)))

cat("\n--- HARNESS ---\n")
cat(sprintf("att=%.10f\natt_se=%.10f\natt_lcl=%.10f\natt_ucl=%.10f\n", att, se, att - z * se, att + z * se))
cat(sprintf("naive_difference=%.10f\n", naive))
cat(sprintf("smd_severity_before=%.10f\nsmd_severity_after=%.10f\n", smd(d), smd(md)))
cat(sprintf("matched_treated=%d\nmatched_pairs=%d\nn=%d\n", sum(md$treated == 1), length(unique(md$subclass)), nrow(d)))
