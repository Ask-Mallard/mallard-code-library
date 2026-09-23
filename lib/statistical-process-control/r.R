# Statistical process control: u-, p- and individuals charts with baseline limits.
#
# 1. PICK THE CHART FOR THE DATA. Counts over varying exposure (infections per catheter-day): u-chart.
#    Proportions over varying denominators (readmissions per discharge): p-chart. One continuous value
#    per period (mean length of stay): individuals (XmR) chart, its spread from the average moving range
#    divided by 1.128. The u- and p-chart limits vary with each month's denominator.
#
# 2. SET LIMITS FROM A BASELINE, THEN APPLY THEM. Limits computed from the first 12 months (phase I) are
#    applied to all 24 (phase II); limits recomputed from data that include a change absorb it. Pinned:
#    3-sigma limits, baseline months 1 to 12.
#
# 3. A SIGNAL IS A SPECIAL CAUSE TO INVESTIGATE, not a significance test. Only the "point beyond a limit"
#    rule is counted here; run rules (for example 8 points on one side of the centre) add sensitivity and
#    false alarms. SPC shows whether a process changed; it does not attribute the change to an
#    intervention, which needs a design (see interrupted-time-series-segmented).

suppressPackageStartupMessages(library(qcc))

d <- read.csv("fixture.csv")
base <- d$month <= 12
later <- !base

u <- qcc(d$infections[base], sizes = d$catheter_days[base], type = "u",
         newdata = d$infections[later], newsizes = d$catheter_days[later], plot = FALSE)
p <- qcc(d$readmissions[base], sizes = d$discharges[base], type = "p",
         newdata = d$readmissions[later], newsizes = d$discharges[later], plot = FALSE)
i <- qcc(d$mean_los[base], type = "xbar.one", std.dev = "MR", newdata = d$mean_los[later], plot = FALSE)

beyond <- function(q) {
  stat <- c(q$statistics, q$newstats)
  lim <- q$limits
  which(stat > lim[, "UCL"] | stat < lim[, "LCL"])
}
first_signal <- function(q) { b <- beyond(q); if (length(b)) b[1] else 0 }
lim_last <- function(q) q$limits[nrow(q$limits), ]

cat(sprintf("u-chart: centre %.3f per 1,000 catheter-days, signals in months %s\n", 1000 * u$center, paste(beyond(u), collapse = ", ")))
cat(sprintf("p-chart: centre %.4f, signals %d; I-chart: centre %.3f, signals %d\n", p$center, length(beyond(p)), i$center, length(beyond(i))))

cat("\n--- HARNESS ---\n")
cat(sprintf("u_centre_per_1000=%.10f\nu_ucl_month24_per_1000=%.10f\n", 1000 * u$center, 1000 * lim_last(u)[["UCL"]]))
cat(sprintf("u_signals=%d\nu_first_signal_month=%d\n", length(beyond(u)), first_signal(u)))
cat(sprintf("p_centre=%.10f\np_ucl_month24=%.10f\np_lcl_month24=%.10f\np_signals=%d\n", p$center, lim_last(p)[["UCL"]], lim_last(p)[["LCL"]], length(beyond(p))))
cat(sprintf("i_centre=%.10f\ni_sigma=%.10f\ni_ucl=%.10f\ni_lcl=%.10f\ni_signals=%d\n", i$center, i$std.dev, lim_last(i)[["UCL"]], lim_last(i)[["LCL"]], length(beyond(i))))
