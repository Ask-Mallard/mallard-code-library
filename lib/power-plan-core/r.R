# Standard power and precision planning, selected by a typed calc key.
#
# The fixture values are examples only. A generated plan must replace them with sourced study
# assumptions. This file deliberately keeps calculation, precision and rule-of-thumb outputs
# distinguishable: a floor is never printed as statistical power.

number <- function(x, default=NA_real_) {
  y <- suppressWarnings(as.numeric(x))
  if (length(y) == 0 || is.na(y)) default else y
}

two_sided_t_power <- function(ncp, df, alpha) {
  critical <- qt(1 - alpha / 2, df)
  pt(-critical, df, ncp=ncp) + pt(critical, df, ncp=ncp, lower.tail=FALSE)
}

exact_t_n <- function(delta, sd, alpha, target, ratio=1, two_groups=FALSE) {
  z <- qnorm(1 - alpha / 2) + qnorm(target)
  if (two_groups) {
    start <- max(2, floor((1 + 1 / ratio) * z^2 * sd^2 / delta^2) - 3)
    for (n1 in start:100000) {
      n2 <- max(2, ceiling(ratio * n1))
      se <- sd * sqrt(1 / n1 + 1 / n2)
      if (two_sided_t_power(abs(delta) / se, n1 + n2 - 2, alpha) >= target)
        return(c(n1=n1, n2=n2, total=n1+n2))
    }
  } else {
    start <- max(2, floor(z^2 * sd^2 / delta^2) - 3)
    for (n in start:100000)
      if (two_sided_t_power(abs(delta) * sqrt(n) / sd, n - 1, alpha) >= target)
        return(c(total=n))
  }
  stop("target power was not reached within the search bound")
}

power_plan <- function(calc, x) {
  alpha <- number(x$alpha, 0.05); target <- number(x$power, 0.8)
  ratio <- number(x$ratio, 1); za <- qnorm(1 - alpha / 2); zb <- qnorm(target)
  if (calc == "twoMeans") return(exact_t_n(number(x$delta), number(x$sd), alpha, target, ratio, TRUE))
  if (calc %in% c("paired", "oneMean")) return(exact_t_n(number(x$delta), number(x$sd), alpha, target))
  if (calc %in% c("twoProportions", "oneProportion")) {
    p1 <- number(x$p1); p2 <- number(x$p2); d <- abs(p1-p2)
    if (calc == "twoProportions") {
      pbar <- (p1 + ratio*p2)/(1+ratio)
      a <- za*sqrt((1+1/ratio)*pbar*(1-pbar))
      b <- zb*sqrt(p1*(1-p1)+p2*(1-p2)/ratio)
      n1 <- (a+b)^2/d^2
      return(c(n1=ceiling(n1),n2=ceiling(ratio*n1),total=ceiling(n1)+ceiling(ratio*n1)))
    }
    n <- (za*sqrt(p1*(1-p1))+zb*sqrt(p2*(1-p2)))^2/d^2
    return(c(total=ceiling(n)))
  }
  if (calc == "correlation") {
    fisher <- atanh(abs(number(x$r)))
    return(c(total=ceiling(((za+zb)/fisher)^2+3)))
  }
  if (calc == "logrank") {
    hr <- number(x$hr); event_rate <- number(x$eventRate); pi1 <- 1/(1+ratio)
    events <- (za+zb)^2/(pi1*(1-pi1)*log(hr)^2)
    return(c(events=ceiling(events),total=ceiling(events/event_rate)))
  }
  if (calc == "logisticEPV") {
    events <- 10*round(number(x$predictors)); return(c(events=events,total=ceiling(events/number(x$eventRate))))
  }
  if (calc == "linearSPV") {
    predictors <- round(number(x$predictors)); return(c(total=max(50+8*predictors,104+predictors)))
  }
  if (calc == "dxAccuracy") {
    se <- number(x$sens); sp <- number(x$spec); prev <- number(x$prev); width <- number(x$width)
    diseased <- ceiling(za^2*se*(1-se)/width^2); non_diseased <- ceiling(za^2*sp*(1-sp)/width^2)
    return(c(diseased=diseased,nonDiseased=non_diseased,total=max(ceiling(diseased/prev),ceiling(non_diseased/(1-prev)))))
  }
  if (calc == "agreementKappa") {
    kappa <- number(x$kappa); prev <- number(x$prev); width <- number(x$width)
    v <- (1-kappa)*((1-kappa)*(1-2*kappa)+kappa*(2-kappa)/(2*prev*(1-prev)))
    return(c(total=max(20,ceiling(za^2*v/width^2))))
  }
  if (calc == "anovaKGroups") {
    groups <- round(number(x$groups)); effect <- number(x$f)
    for (per_group in 2:100000) {
      total <- groups*per_group; critical <- qf(1-alpha,groups-1,total-groups)
      attained <- pf(critical,groups-1,total-groups,ncp=effect^2*total,lower.tail=FALSE)
      if (attained >= target) return(c(perGroup=per_group,total=total))
    }
  }
  if (calc == "proportionPrecision") {
    p <- number(x$p); width <- number(x$width); icc <- number(x$icc,0)
    clusters <- number(x$clusters); cluster_size <- number(x$clusterSize)
    effective <- za^2*p*(1-p)/width^2
    if (!is.na(cluster_size)) exact <- effective*(1+(cluster_size-1)*icc)
    else if (!is.na(clusters)) {
      denominator <- 1-effective*icc/clusters
      if (denominator <= 0) stop("more clusters are required to reach this precision")
      exact <- effective*(1-icc)/denominator
    } else exact <- effective
    total <- ceiling(exact-1e-9)
    if (!is.na(clusters)) total <- max(total,clusters)
    return(c(total=total))
  }
  if (calc == "noninferiorityProportions") {
    p_ref <- number(x$pRef); p_new <- number(x$pNew,p_ref); margin <- number(x$margin)
    one_alpha <- number(x$alpha,0.025); ratio <- number(x$ratio,1)
    expected <- if (tolower(x$worseIs) == "lower") p_ref-p_new else p_new-p_ref
    room <- margin-expected
    n_ref <- ceiling((qnorm(1-one_alpha)+qnorm(target))^2*(p_ref*(1-p_ref)+p_new*(1-p_new)/ratio)/room^2)
    n_new <- ceiling(n_ref*ratio)
    return(c(nReference=n_ref,nNew=n_new,total=n_ref+n_new))
  }
  stop(paste("unsupported calc", calc))
}

d <- read.csv("fixture.csv", stringsAsFactors=FALSE)
results <- list()
for (calc in unique(d$calc)) {
  rows <- d[d$calc == calc,]
  inputs <- as.list(setNames(rows$value, rows$key))
  values <- power_plan(calc, inputs)
  for (name in names(values)) results[[paste(calc,name,sep="_")]] <- values[[name]]
}

cat("These are analyzable sample sizes before attrition or missing-data inflation.\n")
cat("logisticEPV and linearSPV are rule-of-thumb floors, not power calculations.\n")
cat("\n--- HARNESS ---\n")
for (name in names(results)) cat(sprintf("%s=%.10f\n",name,results[[name]]))
