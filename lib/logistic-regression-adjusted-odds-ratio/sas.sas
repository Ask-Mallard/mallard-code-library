/* Logistic regression for an ADJUSTED ODDS RATIO, with profile-likelihood and Wald intervals

   NOT EXECUTED. No SAS engine is licensed for this repository, so unlike the R and Python files
   beside it nothing has ever run this one. It is a reviewed reference implementation and it
   carries neither the agreement claim nor the recovery claim. Read it as a starting point.

   THE PINNED DEFAULT HERE IS THE MOST DANGEROUS ONE IN THE ENTRY. PROC LOGISTIC models the LOWER
   ordered response level. With a 0/1 outcome and no event= option it therefore models P(outcome=0)
   and every odds ratio in the listing is the reciprocal of the one intended -- a 1.87 prints as
   0.53, the output looks entirely normal, and nothing in it says which level was modelled. */

proc import datafile="fixture.csv" out=fixture dbms=csv replace;
    getnames=yes;
run;

/* event='1' IS LOAD-BEARING. See above. The repository's structural lint requires it on any PROC
   LOGISTIC model statement precisely because its absence is invisible in the results. */
proc logistic data=fixture;
    model outcome(event='1') = exposed covariate;

    /* PROFILE-LIKELIHOOD intervals, matching R's confint() and the Python file's brentq root.
       SAS's DEFAULT IS WALD: without clodds=pl the listing reports Wald limits under the heading
       "Odds Ratio Estimates", which is a different interval from the one the R file prints under
       the same name. Both are requested here so neither is implied. */
    oddsratio exposed / cl=pl;
    ods output OddsRatiosWald=wald_limits;
run;

/* The Wald limits for comparison, from the same fit. */
proc print data=wald_limits;
run;
