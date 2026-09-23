# Decision curve analysis: the net benefit of acting on a prediction model's risk.
#
# 1. THE QUESTION IS WHETHER USING THE MODEL IMPROVES DECISIONS, not whether it discriminates. At a
#    threshold risk t (treat when predicted risk exceeds t), net benefit = TP/n - FP/n x t / (1 - t): true
#    positives, less false positives weighted by the odds at which a clinician is indifferent
#    (Vickers and Elkin 2006).
#
# 2. COMPARE WITH THE DEFAULT STRATEGIES: treat everyone (net benefit = prevalence - (1 - prevalence) x
#    t / (1 - t)) and treat no one (0). A model is useful at a threshold only if it beats both.
#
# 3. CHOOSE THE RANGE OF THRESHOLDS BEFOREHAND, from the clinical decision (here 10% to 30%), and plot
#    the whole curve over it. Net benefit uses the predicted risks as given, so a miscalibrated model is
#    penalized; that is intended. Intervals, if wanted, come from a bootstrap. The dcurves package draws
#    the same quantity.

d <- read.csv("fixture.csv")
n <- nrow(d)
prev <- mean(d$event)
nb <- function(t) {
  treat <- d$predicted_risk > t
  sum(treat & d$event == 1) / n - sum(treat & d$event == 0) / n * t / (1 - t)
}
nb_all <- function(t) prev - (1 - prev) * t / (1 - t)

for (t in c(0.1, 0.2, 0.3)) cat(sprintf("threshold %.1f: model %.4f, treat all %.4f, treat none 0\n", t, nb(t), nb_all(t)))

cat("\n--- HARNESS ---\n")
for (t in c(0.1, 0.2, 0.3)) {
  key <- sub(".", "", sprintf("%.1f", t), fixed = TRUE)
  cat(sprintf("net_benefit_model_t%s=%.10f\nnet_benefit_all_t%s=%.10f\n", key, nb(t), key, nb_all(t)))
}
cat(sprintf("prevalence=%.10f\nn=%d\n", prev, n))
