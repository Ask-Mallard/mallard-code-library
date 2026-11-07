# Cohen's and quadratic-weighted kappa for two raters

**R calculator reference: PENDING.** Python recovery passed in this task. A saved R 4.6.1
reference and R/Python agreement verification remain pending installation of `irr 0.85`
from the Posit 2026-09-01 snapshot. The build environment predates access to the package
download host. A subsequent task in the updated environment must produce the R values
and run both languages through `harness/check.py`. No R reference values are supplied here;
the existing recovery targets below describe the synthetic data-generating process.

Two radiologists grade 300 patients on a 1 to 4 scale.

| This fixture | Truth | Estimate |
|---|---|---|
| percent agreement | | 0.51 |
| **Cohen's kappa** | 0.39 | 0.33 (SE 0.040) |
| **quadratic-weighted kappa** | 0.74 | 0.67 (SE 0.030) |

- **Percent agreement is not agreement beyond chance.** Kappa = (observed − chance) / (1 − chance).
- **An ordinal scale needs weights.** Unweighted kappa counts a 1-versus-2 disagreement the same as
  1-versus-4; quadratic weights, 1 − (i − j)² / (k − 1)², give partial credit. Name the weights.
- **The SE** for an interval is the large-sample one of Fleiss, Cohen and Everitt (1969); `irr::kappa2`
  reports only a test of kappa = 0.
- Kappa depends on how common each category is; report the marginal distributions beside it.

## Verification

| Engine | Status |
|---|---|
| R (irr 0.85 for the kappas; SE written out) | executed in CI |
| Python (written out) | executed in CI |
