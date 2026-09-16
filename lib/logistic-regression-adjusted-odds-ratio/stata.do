* Logistic regression for an ADJUSTED ODDS RATIO
*
* NOT EXECUTED. No Stata licence is held for this repository, so nothing has run this file. It
* carries neither the agreement claim nor the recovery claim.
*
* THE PINNED DEFAULT: `logit` reports COEFFICIENTS (log odds ratios) and `logistic` reports ODDS
* RATIOS, from the identical fit. Reading a `logit` table as though it held odds ratios understates
* every effect -- 0.63 read as an odds ratio rather than exp(0.63) = 1.87. Both commands appear
* below so the relationship is on the page rather than in the reader's memory.

import delimited "fixture.csv", clear

* Coefficients: the log odds ratio, directly comparable with the R and Python harness output.
logit outcome exposed covariate

* The same fit reported as odds ratios. `or` is what changes the scale of the table.
logit outcome exposed covariate, or

* PROFILE-LIKELIHOOD intervals are NOT what the table above prints: Stata reports Wald limits, as
* its default and with no option to change them here. The likelihood-ratio interval must be built
* from the profile explicitly, and official Stata has no built-in command for it -- the closest
* supported route is a likelihood-ratio test at candidate values, or the user-written `pllf`.
* This file does not claim a profile interval it cannot produce. The R and Python files report one.
