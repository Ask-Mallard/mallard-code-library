# Cluster-level summary analysis for a cluster trial with few clusters
#
# Summarize each cluster (its infection proportion), then compare the summaries between arms with a
# t-test on (number of clusters - 2) degrees of freedom. With few clusters (here 24) this is valid where
# individual-level sandwich or GEE standard errors are too small. It estimates the difference in the
# AVERAGE CLUSTER proportion, which equals the patient-level difference when clusters are equal in size.

d <- read.csv("fixture.csv")
stopifnot(!anyNA(d), all(d$infection %in% c(0, 1)))

# PINNED: one row per CLUSTER, and var.equal = TRUE, the classical two-sample t on k - 2 df. The unit of
# analysis is the cluster; a test on the individual rows treats 1440 correlated patients as independent.
clusters <- aggregate(infection ~ cluster + arm, data = d, FUN = mean)
fit <- t.test(infection ~ factor(arm, levels = c(1, 0)), data = clusters, var.equal = TRUE, conf.level = 0.95)
diff <- unname(fit$estimate[1] - fit$estimate[2])

cat(sprintf("%d clusters; difference in mean cluster proportion (intervention - control) %.4f, 95%% %.4f to %.4f, %d df\n",
            nrow(clusters), diff, fit$conf.int[1], fit$conf.int[2], as.integer(fit$parameter)))

cat("\n--- HARNESS ---\n")
cat(sprintf("difference=%.10f\ndifference_lcl=%.10f\ndifference_ucl=%.10f\n", diff, fit$conf.int[1], fit$conf.int[2]))
cat(sprintf("df=%d\nclusters=%d\nn=%d\n", as.integer(fit$parameter), nrow(clusters), nrow(d)))
