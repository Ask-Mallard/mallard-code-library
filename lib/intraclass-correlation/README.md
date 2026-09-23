# Intraclass correlation for inter-rater reliability, ICC(A,1)

"The ICC" is six coefficients (Shrout and Fleiss 1979; McGraw and Wong 1996). Name three choices:

- **model**: one-way (each subject rated by different raters) or two-way (the same raters rate every
  subject; random if they stand for raters in general);
- **type**: **absolute agreement** (a rater who reads 5 mm larger counts as disagreeing) or consistency
  (only the ranking matters);
- **unit**: a single rater, or the average of the k raters.

This entry pins **two-way random, absolute agreement, single rater: ICC(A,1)**, the usual answer to "can
one clinician's measurement be trusted", with McGraw and Wong's interval, and reports ICC(C,1) beside it.

| This fixture (60 subjects, 4 raters) | Truth | Estimate |
|---|---|---|
| **ICC(A,1)** | 0.719 | 0.622 (95% CI 0.383 to 0.774) |
| ICC(C,1) | 0.800 | 0.743 |

The ICC rises with the spread of the subjects for the same measurement error; report the between-subject
SD, and for comparing two methods use Bland-Altman limits (`bland-altman-limits`).

## Verification

| Engine | Status |
|---|---|
| R (irr 0.85) | executed in CI |
| Python (irr's formulas written out) | executed in CI |
