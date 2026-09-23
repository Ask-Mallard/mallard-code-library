# Pearson and Spearman correlation with intervals

The correlation between two continuous measurements, with its interval.

| Coefficient | Measures | Interval |
|---|---|---|
| Pearson r | linear association | Fisher z, SE 1/√(n − 3) |
| Spearman ρ | monotonic association, on ranks; resists outliers | Fisher z, SE 1.06/√(n − 3) (Fieller, Hartley and Pearson 1957) |

## The default this entry pins

Neither R's `cor.test` nor scipy's `spearmanr` gives an interval for Spearman's ρ. Both files build it
on Fisher's z with the Fieller standard error; using Pearson's 1/√(n − 3) instead gives a narrower
interval, which the controls show is a different number.

## What a correlation is not

- **Not agreement.** Two methods measuring the same thing can correlate perfectly while one reads 10
  units high. Agreement needs Bland-Altman limits of agreement or an intraclass correlation.
- **Not causation**, and not the slope: r says how tightly points follow a line, not how steep it is.

## The fixture

150 pairs from a bivariate normal with ρ = 0.5 (BMI and systolic pressure). For a bivariate normal,
Spearman's value is (6/π) asin(ρ/2) = 0.4826. Observed: r 0.449, ρ 0.414.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
