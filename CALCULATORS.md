# Calculator formulas and library coverage

All 17 calculator categories listed below already have runnable entries on `main`.
There are no missing entries to add. Pearson and Spearman share one entry; diagnostic
accuracy and predictive values use two complementary entries. Each linked entry has
its README, metadata, reproducible synthetic fixture, R and Python code, and expected results.
The formulas and examples are covered by this repository's [Apache-2.0 license](LICENSE).

| Calculator | Existing entry | Method and standard reference |
|---|---|---|
| Sample size and power | [power-plan-core](lib/power-plan-core/) | Noncentral t and F power, normal approximations for proportions, Fisher-z correlation planning, event-based log-rank planning, and explicitly labeled precision or model-size floors. Cohen (1988); Schoenfeld (1981); Chow et al. (2017); Fleiss et al. (2003); Buderer (1996). Model-size floors: Green (1991); Peduzzi et al. (1996). |
| Descriptive summary: mean and SD, median and IQR | [descriptive-summary-table-one](lib/descriptive-summary-table-one/) | Arithmetic mean; sample SD with divisor n − 1; median and quartiles using R quantile type 7. IQR width is Q3 − Q1. Hyndman and Fan (1996). |
| Chi-square 2x2 with Yates | [chi-square-two-by-two](lib/chi-square-two-by-two/) | Sum of squared observed-minus-expected counts divided by expected counts, with the continuity adjustment capped at the observed deviation. The entry reports both corrected and uncorrected statistics; `chisq_yates` is the corrected calculator result. Yates (1934). |
| Fisher's exact 2x2 | [fisher-exact-conditional-odds-ratio](lib/fisher-exact-conditional-odds-ratio/) | Conditional hypergeometric two-sided test; conditional maximum-likelihood odds ratio and exact conditional interval, distinguished from ad/bc. Fisher (1935). |
| Independent t-test (Welch) | [welch-two-sample-t-test](lib/welch-two-sample-t-test/) | Mean difference divided by sqrt(s1²/n1 + s0²/n0), with Welch-Satterthwaite degrees of freedom and a t interval. Welch (1947). |
| Mann-Whitney U | [mann-whitney-hodges-lehmann](lib/mann-whitney-hodges-lehmann/) | Exact rank-sum test on the tie-free fixture; median of all between-group differences and an interval from inversion of U. Mann and Whitney (1947); Hodges and Lehmann (1963). |
| Paired t-test | [paired-t-test](lib/paired-t-test/) | One-sample t inference on after-minus-before differences: mean change ± t times SD(change)/sqrt(n). Student (1908). |
| Wilcoxon signed-rank | [wilcoxon-signed-rank-hodges-lehmann](lib/wilcoxon-signed-rank-hodges-lehmann/) | Exact signed-rank test; median of Walsh averages and an interval from test inversion. Symmetric changes, with no zeros or tied absolute changes in the fixture. Wilcoxon (1945); Hodges and Lehmann (1963). |
| Pearson correlation | [correlation-pearson-spearman](lib/correlation-pearson-spearman/) | Standardized covariance, with a Fisher-z interval using SE 1/sqrt(n − 3). Pearson (1895); Fisher (1921). |
| Spearman correlation | [correlation-pearson-spearman](lib/correlation-pearson-spearman/) | Pearson correlation of ranks; approximate Fisher-z interval using SE 1.06/sqrt(n − 3). Spearman (1904); Fieller, Hartley and Pearson (1957). |
| Diagnostic test 2x2 | [diagnostic-accuracy-wilson](lib/diagnostic-accuracy-wilson/), [predictive-values-likelihood-ratios](lib/predictive-values-likelihood-ratios/) | Sensitivity TP/(TP+FN), specificity TN/(TN+FP), PPV TP/(TP+FP), NPV TN/(TN+FN), each with a Wilson interval. Likelihood ratios use log-scale intervals. Wilson (1927); Simel, Samsa and Matchar (1991). |
| One-way ANOVA | [one-way-anova-tukey](lib/one-way-anova-tukey/) | Between-group divided by within-group mean square; categorical groups and equal-variance F test, followed by Tukey intervals. Fisher (1925); Tukey (1949). |
| Kruskal-Wallis | [kruskal-wallis](lib/kruskal-wallis/) | Rank-based H statistic, tie correction and chi-square approximation; group medians reported separately. Kruskal and Wallis (1952). |
| Multiple-comparison adjustment | [multiplicity-adjustment](lib/multiplicity-adjustment/) | Holm step-down and Hochberg step-up family-wise error control; Benjamini-Hochberg false discovery rate control. Holm (1979); Hochberg (1988); Benjamini and Hochberg (1995). |
| Kaplan-Meier | [kaplan-meier-logrank](lib/kaplan-meier-logrank/) | Product of (1 − events/at risk), with censoring and ties; Greenwood variance on the log-log scale; log-rank comparison. Kaplan and Meier (1958); Kalbfleisch and Prentice (2002). |
| Cohen's kappa | [weighted-kappa-agreement](lib/weighted-kappa-agreement/) | (Observed agreement − chance agreement)/(1 − chance agreement), with unweighted and quadratic-weighted results and the large-sample Fleiss-Cohen-Everitt SE. Cohen (1960); Fleiss, Cohen and Everitt (1969). |
| Bland-Altman | [bland-altman-limits](lib/bland-altman-limits/) | Mean difference ± normal 0.975 quantile times SD(differences); t interval for bias and approximate t intervals for limits with SE sqrt(3s²/n). Bland and Altman (1986, 1999). |

## Numerical references

Fourteen entries have saved R 4.6.1 references and passed the R/Python agreement check.
Three saved references are **pending**:

| Entry | Required package | Current verification |
|---|---|---|
| `diagnostic-accuracy-wilson` | `binom 1.1-2` | Python recovery passed; R reference and R/Python agreement pending |
| `kaplan-meier-logrank` | `survival 3.8-11` | Python recovery passed; R reference and R/Python agreement pending |
| `weighted-kappa-agreement` | `irr 0.85` | Python recovery passed; R reference and R/Python agreement pending |

Their READMEs and `expected.json` files explicitly record this status. The build environment
predates access to Posit's package download host, `rspm-sync.rstudio.com`. A subsequent task
in the updated environment must install the declared versions and record the actual R outputs.
Existing recovery targets, tolerances, execution declarations and statistical assertions are retained.
The reference test remains strict: it reports these three pending references as failures, while
the other fourteen pass. No pending result is counted as tested agreement.

The `r_reference` section of an entry's `expected.json`, when present, records R 4.6.1
output from its committed fixture, the fixture's SHA-256 digest, and declared R package
versions. These fixed numerical references supplement the existing simulated-parameter
recovery checks. A missing section means the R reference has not been recorded.

Run `python harness/calculator_examples.test.py` to execute both languages, run
`harness/check.py`, and compare every emitted numeric key with its saved R reference.
The check requires references for all entries in `harness/calculator_coverage.json`;
missing references or missing engines fail. Numerical comparisons use the existing
agreement tolerances, including the stated root-finding allowance for Fisher's exact interval.
The test is also discovered by `harness/new_examples.test.py`.

Wilson and Clopper-Pearson intervals apply to a single binomial proportion, including
each sensitivity, specificity or predictive value with its own denominator. They do not
apply to a mean difference, odds ratio, rank-based shift, correlation, survival curve,
kappa, or limit of agreement. The table states the interval appropriate to each estimator.
Power and precision planning use the assumptions documented by their entry; a rule-of-thumb
model-size floor does not establish power.

## R installation and verification

R 4.6.1 was built natively on Debian 13 from the official
`https://cloud.r-project.org/src/base/R-4/R-4.6.1.tar.gz` release archive
(SHA-256 `4da6e61d2c0aac5f14a2e7e432cb5fcc269efe83da4293050ba7f03dff4e2cf4`).
The build used Debian GCC/GFortran 14.2.0-19 and native development headers in a
user-owned directory, with `--enable-R-shlib --with-x=no --with-blas=no --with-lapack=no`.
R uses its own BLAS and LAPACK. The Fortran runtime libraries were linked explicitly;
the installed interpreter's shared-library dependencies resolve on Debian 13.
Version and numerical smoke checks and R's full `make check` passed.

The Posit 2026-09-01 package index matches all 23 R package versions declared by the
library. Installing those packages remains pending; bundled recommended packages do
not satisfy later snapshot pins. In the updated environment, use source packages from
`https://packagemanager.posit.co/cran/2026-09-01`, the same dated CRAN snapshot as CI's
Ubuntu noble endpoint, and run `python harness/metadata.py --installed` to verify every pin.
No interpreter version, package pin, existing tolerance or assertion was changed to
make setup pass.

All 17 calculator Python implementations passed recovery. Fourteen entries also passed
both-language `harness/check.py` comparisons and comparisons with their saved R outputs.
Metadata, comparator, lint and power-plan controls passed, as did structural lint over
all 89 entries and 15 independent example-control suites, including the base-R Bayesian
controls. The full `new_examples.test.py` run reaches the three pending reference failures.
The installed-package check and IPTW scope check remain blocked by unavailable R packages.

## References

- Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *JRSS B* 57:289–300.
- Bland JM, Altman DG (1986). Statistical methods for assessing agreement between two methods of clinical measurement. *Lancet* 1:307–310. (1999). Measuring agreement in method comparison studies. *Statistical Methods in Medical Research* 8:135–160.
- Buderer NMF (1996). Statistical methodology: I. Incorporating the prevalence of disease into the sample size calculation for sensitivity and specificity. *Academic Emergency Medicine* 3:895–900.
- Chow SC, Shao J, Wang H, Lokhnygina Y (2017). *Sample Size Calculations in Clinical Research*, 3rd ed. CRC Press.
- Cohen J (1960). A coefficient of agreement for nominal scales. *Educational and Psychological Measurement* 20:37–46. (1988). *Statistical Power Analysis for the Behavioral Sciences*, 2nd ed. Lawrence Erlbaum.
- Fieller EC, Hartley HO, Pearson ES (1957). Tests for rank correlation coefficients. I. *Biometrika* 44:470–481.
- Fisher RA (1921). On the probable error of a coefficient of correlation deduced from a small sample. *Metron* 1:3–32. (1925). *Statistical Methods for Research Workers*. Oliver and Boyd. (1935). *The Design of Experiments*. Oliver and Boyd.
- Fleiss JL, Cohen J, Everitt BS (1969). Large sample standard errors of kappa and weighted kappa. *Psychological Bulletin* 72:323–327.
- Fleiss JL, Levin B, Paik MC (2003). *Statistical Methods for Rates and Proportions*, 3rd ed. Wiley.
- Green SB (1991). How many subjects does it take to do a regression analysis? *Multivariate Behavioral Research* 26:499–510.
- Hochberg Y (1988). A sharper Bonferroni procedure for multiple tests of significance. *Biometrika* 75:800–802.
- Hodges JL, Lehmann EL (1963). Estimates of location based on rank tests. *Annals of Mathematical Statistics* 34:598–611.
- Holm S (1979). A simple sequentially rejective multiple test procedure. *Scandinavian Journal of Statistics* 6:65–70.
- Hyndman RJ, Fan Y (1996). Sample quantiles in statistical packages. *The American Statistician* 50:361–365.
- Kalbfleisch JD, Prentice RL (2002). *The Statistical Analysis of Failure Time Data*, 2nd ed. Wiley.
- Kaplan EL, Meier P (1958). Nonparametric estimation from incomplete observations. *JASA* 53:457–481.
- Kruskal WH, Wallis WA (1952). Use of ranks in one-criterion variance analysis. *JASA* 47:583–621.
- Mann HB, Whitney DR (1947). On a test of whether one of two random variables is stochastically larger than the other. *Annals of Mathematical Statistics* 18:50–60.
- Pearson K (1895). Notes on regression and inheritance in the case of two parents. *Proceedings of the Royal Society of London* 58:240–242.
- Peduzzi P, Concato J, Kemper E, Holford TR, Feinstein AR (1996). A simulation study of the number of events per variable in logistic regression analysis. *Journal of Clinical Epidemiology* 49:1373–1379.
- Schoenfeld D (1981). The asymptotic properties of nonparametric tests for comparing survival distributions. *Biometrika* 68:316–319.
- Simel DL, Samsa GP, Matchar DB (1991). Likelihood ratios with confidence: sample size estimation for diagnostic test studies. *Journal of Clinical Epidemiology* 44:763–770.
- Spearman C (1904). The proof and measurement of association between two things. *The American Journal of Psychology* 15:72–101.
- Student (1908). The probable error of a mean. *Biometrika* 6:1–25.
- Tukey JW (1949). Comparing individual means in the analysis of variance. *Biometrics* 5:99–114.
- Welch BL (1947). The generalization of Student's problem when several different population variances are involved. *Biometrika* 34:28–35.
- Wilcoxon F (1945). Individual comparisons by ranking methods. *Biometrics Bulletin* 1:80–83.
- Wilson EB (1927). Probable inference, the law of succession, and statistical inference. *JASA* 22:209–212.
- Yates F (1934). Contingency tables involving small numbers and the chi-square test. *Supplement to the JRSS* 1:217–235.
