# Inter-rater reliability: the intraclass correlation, with the model, type and unit named.
#
# 1. "THE ICC" IS SIX DIFFERENT COEFFICIENTS (Shrout and Fleiss 1979; McGraw and Wong 1996). Name three
#    choices, because they change the number:
#    - model: one-way (each subject rated by different raters) or two-way (the same raters rate every
#      subject; raters random if they stand for raters in general);
#    - type: ABSOLUTE AGREEMENT (a rater who reads 5 mm larger than the others counts as disagreement) or
#      CONSISTENCY (only the ranking matters);
#    - unit: a single rater's measurement, or the average of the k raters.
#    Pinned here: two-way random, absolute agreement, single rater, ICC(A,1), the usual choice for "can a
#    single clinician's measurement be trusted". Consistency is reported beside it.
#
# 2. THE INTERVAL for ICC(A,1) is McGraw and Wong's F-based one with a Satterthwaite-type df, as irr::icc
#    computes it.
#
# 3. The ICC depends on how much subjects vary: the same measurement error gives a higher ICC in a more
#    heterogeneous sample. Report the between-subject SD beside it, and for method comparison use
#    Bland-Altman limits (bland-altman-limits).

suppressPackageStartupMessages(library(irr))

d <- read.csv("fixture.csv")
ratings <- as.matrix(d[, grep("^rater_", names(d))])
a <- icc(ratings, model = "twoway", type = "agreement", unit = "single")
cc <- icc(ratings, model = "twoway", type = "consistency", unit = "single")

cat(sprintf("ICC(A,1) %.3f (95%% CI %.3f to %.3f); ICC(C,1) %.3f\n", a$value, a$lbound, a$ubound, cc$value))

cat("\n--- HARNESS ---\n")
cat(sprintf("icc_agreement=%.10f\nicc_agreement_lcl=%.10f\nicc_agreement_ucl=%.10f\n", a$value, a$lbound, a$ubound))
cat(sprintf("icc_consistency=%.10f\nicc_consistency_lcl=%.10f\nicc_consistency_ucl=%.10f\n", cc$value, cc$lbound, cc$ubound))
cat(sprintf("subjects=%d\nraters=%d\n", nrow(ratings), ncol(ratings)))
