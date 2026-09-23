# Target trial emulation with a grace period: clone, censor, weight.
#
# 1. SPECIFY THE TRIAL FIRST. Eligibility and time zero (here: eligible at interval 0), the strategies
#    ("start within 3 intervals" against "never start"), the outcome (an event within 10 intervals), the
#    contrast (per-protocol risk difference). Time zero is the same for both strategies, which is what
#    removes immortal time.
#
# 2. CLONE. At time zero a patient's data are consistent with both strategies (nobody has to start on
#    day one), so each patient is copied into both arms.
#
# 3. CENSOR each clone when its data stop agreeing with its strategy: the "never" clone when the patient
#    starts; the "grace" clone at the end of the grace period if the patient has not started by then.
#
# 4. WEIGHT, because that censoring depends on severity, which also drives the outcome. Each uncensored
#    clone is weighted by the inverse probability of remaining uncensored given its history, from a
#    pooled logistic model of starting. Here: the never arm by prod 1 / (1 - P(start)) over the grace
#    intervals, and grace-arm patients who start in the last grace interval by 1 / P(start).
#
# 5. ESTIMATE the risk in each arm by the weighted discrete-time Kaplan-Meier over the 10 intervals.
#    Confidence intervals need a bootstrap of the whole procedure (not run here: the entry is checked
#    against its truth only).

GRACE <- 3
INTERVALS <- 10

d <- read.csv("fixture.csv")
d <- d[order(d$id, d$interval), ]
start_at <- tapply(ifelse(d$treated == 1, d$interval, Inf), d$id, min)
d$s <- start_at[as.character(d$id)]

# The start (censoring) model: pooled logistic among patient-intervals still eligible to start.
elig <- d$interval < GRACE & d$s >= d$interval
d$start <- as.integer(d$s == d$interval)
ps <- glm(start ~ severity, data = d[elig, ], family = binomial)
d$p <- NA_real_
d$p[elig] <- fitted(ps)

# Never arm: rows before the patient starts; weight = cumulative 1 / (1 - P(start)) over grace intervals.
nv <- d[d$interval < d$s, ]
nv$w <- ave(ifelse(nv$interval < GRACE, 1 / (1 - nv$p), 1), nv$id, FUN = cumprod)

# Grace arm: patients who started in intervals 0 or 1 (weight 1), those who started in the last grace
# interval (weight 1 / P(start) from then on), and non-starters until the end of the grace period.
gr <- d[d$s < GRACE - 1 | d$s == GRACE - 1 | d$interval < GRACE - 1, ]
p_last <- with(d[d$interval == GRACE - 1 & d$s == GRACE - 1, ], setNames(p, id))
gr$w <- ifelse(gr$s == GRACE - 1 & gr$interval >= GRACE - 1, 1 / p_last[as.character(gr$id)], 1)

risk <- function(x, w) {
  h <- sapply(0:(INTERVALS - 1), function(k) {
    at <- x$interval == k
    if (!any(at)) return(0)
    sum(w[at] * x$event[at]) / sum(w[at])
  })
  1 - prod(1 - h)
}

r_grace <- risk(gr, gr$w)
r_never <- risk(nv, nv$w)
rd_unweighted <- risk(gr, rep(1, nrow(gr))) - risk(nv, rep(1, nrow(nv)))

# The immortal-time comparison: patients who ever started against those who never did.
ever <- tapply(d$treated, d$id, max)
any_event <- tapply(d$event, d$id, max)
rd_ever_never <- mean(any_event[ever == 1]) - mean(any_event[ever == 0])

cat(sprintf("10-interval risk: start within grace %.3f, never start %.3f, difference %.3f\n",
            r_grace, r_never, r_grace - r_never))

cat("\n--- HARNESS ---\n")
cat(sprintf("risk_grace=%.10f\nrisk_never=%.10f\nrisk_difference=%.10f\n", r_grace, r_never, r_grace - r_never))
cat(sprintf("risk_difference_unweighted=%.10f\nrisk_difference_ever_vs_never=%.10f\n", rd_unweighted, rd_ever_never))
cat(sprintf("max_weight=%.10f\npatients=%d\n", max(c(gr$w, nv$w)), length(unique(d$id))))
