# Network meta-analysis (frequentist, random effects)

Trials of placebo (A) against B and against C also say something about C against B, through A. A
network meta-analysis combines that indirect evidence with the direct B-versus-C trials.

| This fixture (15 two-arm trials) | Truth (log OR) | Estimate (random effects) |
|---|---|---|
| B vs A | −0.50 | −0.38 (SE 0.11) |
| C vs A | −0.80 | −0.78 (SE 0.13) |
| C vs B | −0.30 | −0.40 (SE 0.14) |
| τ², Q, between-design Q | | 0.044, 21.9, 0.95 |

For C vs B the network's common-effect SE is 0.105, against 0.150 from the four direct trials alone.

## What a network needs

- **Consistency**: direct and indirect evidence agree. Report the between-design (inconsistency) Q.
- **Transitivity**: the trials are similar in the things that modify the effect. A clinical judgement.
- Rankings (P-scores, SUCRA) come after effects with intervals, not instead of them.

## The model

`netmeta`'s model (Rücker 2012): for two-arm trials, weighted least squares on the trial contrasts; a
common τ² by the generalised DerSimonian-Laird estimator; random effects. Both files write it out and
reproduce `netmeta` 3.7.0 to 1e-10 on this fixture; `netmeta` itself cannot load on the CI runner (its
dependency igraph needs a missing system library). With multi-arm trials, use `netmeta`.

## Verification

| Engine | Status |
|---|---|
| R (written out; matches netmeta locally) | executed in CI |
| Python (written out) | executed in CI |
