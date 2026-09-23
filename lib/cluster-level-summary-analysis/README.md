# Cluster-level summary analysis for few clusters

A cluster-randomized trial with few clusters (under about 30). Individual-level sandwich and GEE standard
errors are too small with so few clusters, and analysing patients as if independent is far worse. The
simple, valid alternative: **summarize each cluster, then compare the cluster summaries** with a t-test
on (clusters − 2) degrees of freedom.

| Analysis | 95% half-width here |
|---|---|
| **cluster-level t-test, 22 df** | **0.068** |
| patients as if independent (1440 rows) | 0.044 (far too narrow) |

## What it estimates

The difference in the **average cluster** proportion. With equal cluster sizes, as here, that equals the
patient-level difference. With unequal sizes, weight by size or use a small-sample-corrected model
(Kenward-Roger or CR2).

## The fixture

24 clinics of 60 patients, 12 per arm; clinic risk = 0.35 (control) or 0.20 (intervention) + a clinic
shift uniform on ±0.12. True difference −0.15; observed −0.171.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
