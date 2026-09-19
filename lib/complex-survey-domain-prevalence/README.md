# Complex survey prevalence and domain estimation

Weighted binary prevalence and Taylor-linearized standard errors for a one-stage stratified PSU design, retaining all sampled PSUs for a domain estimate. This is a limited worked example, not a general survey-analysis engine.

Declare strata, PSUs and inverse-inclusion-probability weights on the full eligible sample before domain restriction.

PSU identifiers repeat between strata; R uses nest=TRUE and Python groups by both identifiers.

With-replacement Taylor variance, no finite-population correction. The example covers one-stage PSU sampling with small sampling fractions, not replicate-weight, calibrated or multistage designs.

Every sampled PSU, including PSUs with zero domain observations, contributes to the stratum variance calculation.

Missing values, nonpositive weights and unresolved lonely PSUs are rejected. Provider-specific interval and reliability rules remain necessary.

R and Python are executed on the same synthetic fixture. Agreement and recovery are separate checks. Generated adaptations and analyses on user data have not been executed by these checks.

Sources: https://r-survey.r-forge.r-project.org/pkgdown/docs/reference/surveysummary.html and https://r-survey.r-forge.r-project.org/survey/example-domain.html
