# Noninferiority analysis of a risk difference against a prespecified margin.
#
# 1. THE CONCLUSION COMES FROM THE CONFIDENCE INTERVAL AND THE MARGIN. The new treatment is noninferior
#    if the lower limit of the two-sided 95% interval for (new - standard) lies above the margin
#    (-0.10 here), equivalently a one-sided test at 2.5%. "No significant difference" in a superiority
#    test is not noninferiority.
#
# 2. THE MARGIN IS SET BEFOREHAND, on clinical and historical grounds (the effect of the standard
#    against placebo, and how much of it may be lost), never chosen after seeing the data.
#
# 3. REPORT BOTH ITT AND PER-PROTOCOL. In a superiority trial intention-to-treat is conservative; in a
#    noninferiority trial it is not, because protocol deviations blur the arms together and pull the
#    difference towards 0, towards a noninferiority conclusion. The conclusion should hold in both.
#
# 4. THE INTERVAL is Newcombe's hybrid score interval (method 10), which behaves well near 0 or 1 and in
#    small samples; the Wald interval is reported beside it. Pinned, written out.

d <- read.csv("fixture.csv")
margin <- -0.10
z <- qnorm(0.975)

wilson <- function(x, n) {
  p <- x / n
  centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(centre - half, centre + half)
}
newcombe <- function(x1, n1, x0, n0) {
  p1 <- x1 / n1; p0 <- x0 / n0
  w1 <- wilson(x1, n1); w0 <- wilson(x0, n0)
  diff <- p1 - p0
  c(diff = diff,
    lcl = diff - sqrt((p1 - w1[1])^2 + (w0[2] - p0)^2),
    ucl = diff + sqrt((w1[2] - p1)^2 + (p0 - w0[1])^2),
    wald_lcl = diff - z * sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0))
}
analyse <- function(s) newcombe(sum(s$success[s$randomized_new == 1]), sum(s$randomized_new == 1),
                                sum(s$success[s$randomized_new == 0]), sum(s$randomized_new == 0))
itt <- analyse(d)
pp <- analyse(d[d$per_protocol == 1, ])

for (nm in c("ITT", "PP")) {
  r <- if (nm == "ITT") itt else pp
  cat(sprintf("%s: difference %.3f (95%% CI %.3f to %.3f); noninferior at margin %.2f: %s\n",
              nm, r[["diff"]], r[["lcl"]], r[["ucl"]], margin, r[["lcl"]] > margin))
}

cat("\n--- HARNESS ---\n")
cat(sprintf("risk_difference_itt=%.10f\nitt_lcl=%.10f\nitt_ucl=%.10f\nitt_wald_lcl=%.10f\n", itt[["diff"]], itt[["lcl"]], itt[["ucl"]], itt[["wald_lcl"]]))
cat(sprintf("risk_difference_pp=%.10f\npp_lcl=%.10f\npp_ucl=%.10f\n", pp[["diff"]], pp[["lcl"]], pp[["ucl"]]))
cat(sprintf("noninferior_itt=%d\nnoninferior_pp=%d\n", as.integer(itt[["lcl"]] > margin), as.integer(pp[["lcl"]] > margin)))
cat(sprintf("n=%d\nn_per_protocol=%d\n", nrow(d), sum(d$per_protocol == 1)))
