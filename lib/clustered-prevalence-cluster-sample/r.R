# A prevalence from a cluster sample, with a design-based 95% interval
#
# Clinics were sampled, then everyone in each sampled clinic. Patients in one clinic resemble each
# other, so the information in the sample is closer to the number of CLINICS than the number of
# patients. A Wilson interval over the patients treats them as independent and is too narrow.
#
# One-stage cluster sample, equal-probability clusters, with-replacement (Taylor-linearized)
# variance, no finite population correction: appropriate when the sampled clusters are a small
# fraction of all clusters. Stratified or multistage designs need their own declaration (see
# complex-survey-domain-prevalence). Complete data only.

library(survey)

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d), all(d$weight > 0), all(d$outcome %in% c(0, 1)))

# PINNED: ids = ~cluster. Omitting it (ids = ~1) declares every patient an independently sampled
# unit, and the standard error comes out too small by the square root of the design effect.
design <- svydesign(ids = ~cluster, weights = ~weight, data = d)
est <- svymean(~outcome, design)

# PINNED: df = degf(design), the number of clusters minus one. confint's default is df = Inf, a
# normal interval, which is too narrow when the clusters are few. With 40 clusters the difference is
# small; with 12 it is not.
ci <- confint(est, df = degf(design))

cat(sprintf("%d patients in %d clusters\n", nrow(d), length(unique(d$cluster))))
cat(sprintf("patient-level prevalence %.4f, design-based SE %.4f, 95%% %.4f to %.4f\n",
            coef(est)[1], SE(est)[1], ci[1], ci[2]))

cat("\n--- HARNESS ---\n")
cat(sprintf("prevalence=%.10f\n", coef(est)[1]))
cat(sprintf("prevalence_se=%.10f\n", SE(est)[1]))
cat(sprintf("prevalence_lcl=%.10f\n", ci[1]))
cat(sprintf("prevalence_ucl=%.10f\n", ci[2]))
cat(sprintf("clusters=%d\n", length(unique(d$cluster))))
cat(sprintf("n=%d\n", nrow(d)))
