# Bland-Altman limits of agreement for method comparison

**Correlation is not agreement.** Here the two methods correlate at 0.96, while the new method reads 2.5
units higher on average and individual differences range over 19 units.

| This fixture (B − A, 150 patients) | Truth | Estimate |
|---|---|---|
| bias | 2.0 | 2.45 (95% CI 1.69 to 3.22) |
| lower limit of agreement | −7.80 | −6.87 (−8.20 to −5.54) |
| upper limit of agreement | 11.80 | 11.77 (10.45 to 13.10) |
| proportional-bias slope | 0 | −0.017 |

- **Limits** = bias ± 1.96 SD of the differences. Whether they are narrow enough is a clinical judgement
  to state in the protocol before the study.
- **Intervals** (Bland and Altman 1999): the bias with a t interval, each limit with SE √(3s²/n). The
  exact method (Carkeet 2015) is slightly wider in small samples.
- **Proportional bias**: regress the differences on the means and plot them; a sloping cloud means the
  limits are not constant across the range.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
