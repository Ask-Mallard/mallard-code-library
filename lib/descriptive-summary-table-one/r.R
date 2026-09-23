# Table 1: baseline characteristics overall and by exposure group, with standardized differences
#
# Each variable is summarized by its TYPE, not its position in a list: mean (SD) for a roughly
# symmetric continuous variable, median [Q1, Q3] for a skewed one, n (%) for binary and categorical.
#
# Standardized mean differences (SMD), not p-values, compare the groups. A p-value in Table 1 tests
# a difference the reader already sees, depends on sample size, and in a randomized trial tests a
# null hypothesis that is true by design (CONSORT advises against it). An absolute SMD above about
# 0.1 is the usual flag for imbalance worth adjusting for.
#
# Complete data. With missing values, report the number missing per variable and state whether
# percentages use the non-missing count.

d <- read.csv("fixture.csv", stringsAsFactors = FALSE)
stopifnot(!anyNA(d), all(d$exposed %in% c(0, 1)))

e1 <- d[d$exposed == 1, ]
e0 <- d[d$exposed == 0, ]

# PINNED: quantile type = 7, R's default, named anyway because statistical packages disagree on the
# definition (SAS and SPSS default to others) and quartiles of a small sample move with it.
q <- quantile(d$los_days, c(0.25, 0.5, 0.75), type = 7, names = FALSE)

# sd() divides by n - 1. The SMD pools the two group SDs as sqrt((s1^2 + s0^2) / 2) (Austin 2009),
# and a binary variable uses p(1 - p) as its variance.
smd_continuous <- function(x1, x0) (mean(x1) - mean(x0)) / sqrt((sd(x1)^2 + sd(x0)^2) / 2)
smd_binary <- function(x1, x0) {
  p1 <- mean(x1); p0 <- mean(x0)
  (p1 - p0) / sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
}

prop <- function(x, level) mean(x == level)

cat(sprintf("N = %d (exposed %d, unexposed %d)\n", nrow(d), nrow(e1), nrow(e0)))
cat(sprintf("Age, mean (SD)          %.1f (%.1f) | %.1f (%.1f) vs %.1f (%.1f), SMD %.2f\n",
            mean(d$age), sd(d$age), mean(e1$age), sd(e1$age), mean(e0$age), sd(e0$age),
            smd_continuous(e1$age, e0$age)))
cat(sprintf("Length of stay, median [IQR] %.2f [%.2f, %.2f]\n", q[2], q[1], q[3]))
cat(sprintf("Diabetes, n (%%)        %d (%.1f) | SMD %.2f\n", sum(d$diabetes),
            100 * mean(d$diabetes), smd_binary(e1$diabetes, e0$diabetes)))
for (lvl in c("never", "former", "current"))
  cat(sprintf("Smoking %-8s n (%%)  %d (%.1f)\n", lvl, sum(d$smoking == lvl), 100 * prop(d$smoking, lvl)))

cat("\n--- HARNESS ---\n")
cat(sprintf("n=%d\nn_exposed=%d\nn_unexposed=%d\n", nrow(d), nrow(e1), nrow(e0)))
cat(sprintf("age_mean=%.10f\nage_sd=%.10f\n", mean(d$age), sd(d$age)))
cat(sprintf("age_mean_exposed=%.10f\nage_mean_unexposed=%.10f\n", mean(e1$age), mean(e0$age)))
cat(sprintf("age_sd_exposed=%.10f\nage_sd_unexposed=%.10f\n", sd(e1$age), sd(e0$age)))
cat(sprintf("age_smd=%.10f\n", smd_continuous(e1$age, e0$age)))
cat(sprintf("los_q1=%.10f\nlos_median=%.10f\nlos_q3=%.10f\n", q[1], q[2], q[3]))
cat(sprintf("diabetes_prop=%.10f\n", mean(d$diabetes)))
cat(sprintf("diabetes_prop_exposed=%.10f\ndiabetes_prop_unexposed=%.10f\n",
            mean(e1$diabetes), mean(e0$diabetes)))
cat(sprintf("diabetes_smd=%.10f\n", smd_binary(e1$diabetes, e0$diabetes)))
cat(sprintf("smoking_never_prop=%.10f\n", prop(d$smoking, "never")))
cat(sprintf("smoking_former_prop=%.10f\n", prop(d$smoking, "former")))
cat(sprintf("smoking_current_prop=%.10f\n", prop(d$smoking, "current")))
