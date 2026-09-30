# Fixed-design power. The ITS is simulated; the stepped wedge uses a validated closed form.

number <- function(x, default=NA_real_) { y <- suppressWarnings(as.numeric(x)); if (length(y)==0||is.na(y)) default else y }
read_inputs <- function(calc) {
  d <- read.csv("fixture.csv",stringsAsFactors=FALSE); d <- d[d$calc==calc,]
  as.list(setNames(d$value,d$key))
}

stepped_wedge_schedule <- function(clusters,periods) {
  steps <- periods-1; base <- clusters%/%steps; extra <- clusters%%steps
  per_step <- base+c(rep(1,extra),rep(0,steps-extra)); rows <- matrix(0,clusters,periods); row <- 0
  for (step in seq_len(steps)) for (i in seq_len(per_step[step])) { row<-row+1; rows[row,(step+1):periods]<-1 }
  list(rows=rows,perStep=per_step)
}

stepped_wedge_variance <- function(clusters,periods,size,icc,cac,sigma2) {
  x <- stepped_wedge_schedule(clusters,periods)$rows
  a <- sigma2*(icc*(1-cac)+(1-icc)/size); c <- sigma2*icc*cac
  vi <- (diag(periods)-c/(a+periods*c)*matrix(1,periods,periods))/a
  s <- colSums(x); info <- sum(vapply(seq_len(clusters),function(i)t(x[i,])%*%vi%*%x[i,],numeric(1)))-as.numeric(t(s)%*%vi%*%s)/clusters
  1/info
}

stepped_wedge_power <- function(x) {
  clusters<-round(number(x$clusters)); periods<-round(number(x$periods)); size<-number(x$clusterPeriodSize)
  icc<-number(x$icc); cac<-number(x$cac,1); alpha<-number(x$alpha,0.05); target<-number(x$power,0.8)
  p0<-number(x$baselineRate); p1<-number(x$targetRate); delta<-p1-p0; p<-(p0+p1)/2; sigma2<-p*(1-p)
  power_at <- function(m) pnorm(abs(delta)/sqrt(stepped_wedge_variance(clusters,periods,m,icc,cac,sigma2))-qnorm(1-alpha/2))
  attained<-power_at(size); ceiling_power<-power_at(1e6); needed<-NA_real_
  if(attained>=target) needed<-ceiling(size) else if(ceiling_power>=target) {
    lo<-ceiling(size); hi<-1e6
    while(hi-lo>1){mid<-floor((lo+hi)/2);if(power_at(mid)>=target)hi<-mid else lo<-mid}
    needed<-hi
  }
  c(power=attained,se=sqrt(stepped_wedge_variance(clusters,periods,size,icc,cac,sigma2)),sizeForTarget=needed,ceiling=ceiling_power)
}

its_power <- function(x) {
  pre<-round(number(x$preTimes)); post_n<-round(number(x$postTimes)); controlled<-round(number(x$series,1))==2
  baseline<-number(x$baselineEvents); control_baseline<-number(x$controlBaselineEvents,baseline)
  level_rr<-number(x$rateRatio); slope_rr<-number(x$slopeRatio,1); theta<-number(x$nbSize,0)
  step_at<-number(x$stepAt); step_rr<-number(x$stepRatio,1); amplitude<-number(x$seasonalAmplitude,0)
  alpha<-number(x$alpha,0.05); reps<-round(number(x$reps,1000)); seed<-round(number(x$seed,20260915))
  time<-seq_len(pre+post_n); is_post<-as.integer(time>pre); tpost<-pmax(0,time-pre)
  step<-if(is.na(step_at))rep(0,length(time))else as.integer(time>=step_at)
  sine<-sin(2*pi*(time-1)/12); cosine<-cos(2*pi*(time-1)/12)
  make_series<-function(A){base<-if(A==1)baseline else control_baseline;mu<-base*exp(log(step_rr)*step+amplitude*cosine+A*(log(level_rr)*is_post+log(slope_rr)*tpost));list(A=A,mu=mu)}
  one<-function(){
    series<-if(controlled)list(make_series(1),make_series(0))else list(make_series(1)); d<-do.call(rbind,lapply(series,function(z)data.frame(y=if(theta>0)rnbinom(length(time),size=theta,mu=z$mu)else rpois(length(time),z$mu),A=z$A,time,is_post,tpost,step,sine,cosine)))
    fit<-try(if(controlled)glm(y~A*time+A:time+step+sine+cosine+A:is_post+A:tpost,family=quasipoisson,data=d)else glm(y~time+step+sine+cosine+is_post+tpost,family=quasipoisson,data=d),silent=TRUE)
    if(inherits(fit,"try-error"))return(NA)
    name<-if(controlled)"A:is_post" else "is_post"; co<-summary(fit)$coefficients
    if(!name%in%rownames(co)||!is.finite(co[name,4]))return(NA)
    co[name,4]<alpha
  }
  set.seed(seed); values<-replicate(reps,one()); successful<-sum(!is.na(values)); rejected<-sum(values,na.rm=TRUE); power<-rejected/successful
  c(power=power,mcse=sqrt(power*(1-power)/successful),attempted=reps,successful=successful,failed=reps-successful,seed=seed)
}

its<-its_power(read_inputs("interruptedTimeSeries")); sw<-stepped_wedge_power(read_inputs("steppedWedge"))
cat(sprintf("ITS attainable power %.3f (Monte Carlo SE %.3f; %d/%d successful fits).\n",its["power"],its["mcse"],its["successful"],its["attempted"]))
cat(sprintf("Stepped-wedge attainable power %.3f; size for target %s.\n",sw["power"],ifelse(is.na(sw["sizeForTarget"]),"not attainable",sw["sizeForTarget"])))
cat("\n--- HARNESS ---\n")
for(name in names(its))cat(sprintf("its_%s=%.10f\n",name,its[[name]]))
for(name in names(sw))cat(sprintf("sw_%s=%.10f\n",name,sw[[name]]))
