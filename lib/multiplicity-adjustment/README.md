# Multiplicity adjustment: Holm, Hochberg and Benjamini-Hochberg

Ten secondary outcomes, three with a real effect. Testing each at 0.05 without adjustment rejects at
least one of the seven null outcomes about 30% of the time (100 simulated trials: 30%).

## Decide what needs adjusting, then choose the error rate

- A single prespecified primary outcome needs no adjustment. A **family** of outcomes, subgroups or arms
  that will each support a claim does; name the family in the protocol.
- **Family-wise error rate** (the chance of any false claim): **Holm** (valid under any dependence;
  uniformly better than Bonferroni) or **Hochberg** (more powerful; needs independence or positive
  dependence). Over 100 simulated trials Holm's family-wise error rate was 6%.
- **False discovery rate** (the expected share of false claims among those made): **Benjamini-Hochberg**,
  for screening. It does not control the family-wise rate (13% in the simulations), by design.

Report adjusted p-values beside the unadjusted ones, with the family and the method named.

## Verification

| Engine | Status |
|---|---|
| R (`p.adjust`) | executed in CI |
| Python (statsmodels `multipletests`) | executed in CI |

All 44 p-values agree exactly; the control file rebuilds Holm and Benjamini-Hochberg by hand.
