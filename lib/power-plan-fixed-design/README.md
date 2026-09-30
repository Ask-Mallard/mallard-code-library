# Fixed-design power

The interrupted-time-series route answers “what power does this fixed time grid have?” by fitting
the prespecified segmented quasi-Poisson analysis inside each seeded negative-binomial replicate.
It reports Monte Carlo uncertainty and failed fits. It does not model residual serial correlation,
so its result is an upper bound when monthly residuals remain autocorrelated.

The stepped-wedge route uses the nested-exchangeable GLS variance for a cross-sectional design with
fixed period effects. It reports attainable power at the fixed clusters, periods and cluster-period
size, plus the cluster-period size needed for the target where one exists. With few clusters, a
small-sample reference distribution can yield less power than the normal approximation shown.
