"""Fixed-design power: seeded ITS simulation and stepped-wedge GLS variance."""

import csv
import math
import numpy as np
import statsmodels.api as sm
from scipy.stats import norm, t


def number(value, default=math.nan):
    try:
        return default if value in (None, "") else float(value)
    except (TypeError, ValueError):
        return default


with open("fixture.csv",newline="",encoding="utf-8") as handle:
    rows=list(csv.DictReader(handle))
inputs={}
for row in rows:
    inputs.setdefault(row["calc"],{})[row["key"]]=row["value"]


def stepped_wedge_schedule(clusters,periods):
    steps=periods-1; base,extra=divmod(clusters,steps)
    per_step=[base+(step<extra) for step in range(steps)]; rows=[]
    for step,count in enumerate(per_step):
        rows.extend([[int(period>step) for period in range(periods)] for _ in range(count)])
    return np.asarray(rows,dtype=float),per_step


def stepped_wedge_variance(clusters,periods,size,icc,cac,sigma2):
    design,_=stepped_wedge_schedule(clusters,periods)
    a=sigma2*(icc*(1-cac)+(1-icc)/size); c=sigma2*icc*cac
    inverse=(np.eye(periods)-c/(a+periods*c)*np.ones((periods,periods)))/a
    sums=design.sum(axis=0)
    information=sum(row@inverse@row for row in design)-sums@inverse@sums/clusters
    return 1/information


def stepped_wedge_power(x):
    clusters=round(number(x["clusters"])); periods=round(number(x["periods"])); size=number(x["clusterPeriodSize"])
    icc=number(x["icc"]); cac=number(x.get("cac"),1); alpha=number(x.get("alpha"),0.05); target=number(x.get("power"),0.8)
    p0=number(x["baselineRate"]); p1=number(x["targetRate"]); delta=p1-p0; p=(p0+p1)/2; sigma2=p*(1-p)
    def power_at(m):
        variance=stepped_wedge_variance(clusters,periods,m,icc,cac,sigma2)
        return norm.cdf(abs(delta)/math.sqrt(variance)-norm.ppf(1-alpha/2))
    attained=power_at(size); ceiling_power=power_at(1_000_000); needed=math.nan
    if attained>=target: needed=math.ceil(size)
    elif ceiling_power>=target:
        low,high=math.ceil(size),1_000_000
        while high-low>1:
            middle=(low+high)//2
            if power_at(middle)>=target: high=middle
            else: low=middle
        needed=high
    variance=stepped_wedge_variance(clusters,periods,size,icc,cac,sigma2)
    return {"power":attained,"se":math.sqrt(variance),"sizeForTarget":needed,"ceiling":ceiling_power}


def its_power(x):
    pre=round(number(x["preTimes"])); post_n=round(number(x["postTimes"])); controlled=round(number(x.get("series"),1))==2
    baseline=number(x["baselineEvents"]); control_baseline=number(x.get("controlBaselineEvents"),baseline)
    level_rr=number(x["rateRatio"]); slope_rr=number(x.get("slopeRatio"),1); theta=number(x.get("nbSize"),0)
    step_at=number(x.get("stepAt")); step_rr=number(x.get("stepRatio"),1); amplitude=number(x.get("seasonalAmplitude"),0)
    alpha=number(x.get("alpha"),0.05); reps=round(number(x.get("reps"),1000)); seed=round(number(x.get("seed"),20260915))
    time=np.arange(1,pre+post_n+1,dtype=float); is_post=(time>pre).astype(float); tpost=np.maximum(0,time-pre)
    step=np.zeros_like(time) if not math.isfinite(step_at) else (time>=step_at).astype(float)
    sine=np.sin(2*np.pi*(time-1)/12); cosine=np.cos(2*np.pi*(time-1)/12)
    rng=np.random.default_rng(seed); rejected=0; failed=0
    for _ in range(reps):
        blocks=[]
        for arm in ([1,0] if controlled else [1]):
            base=baseline if arm else control_baseline
            mu=base*np.exp(np.log(step_rr)*step+amplitude*cosine+arm*(np.log(level_rr)*is_post+np.log(slope_rr)*tpost))
            if theta>0:
                gamma=rng.gamma(shape=theta,scale=mu/theta); outcome=rng.poisson(gamma)
            else: outcome=rng.poisson(mu)
            if controlled:
                design=np.column_stack([np.ones(len(time)),np.full(len(time),arm),time,arm*time,step,sine,cosine,arm*is_post,arm*tpost])
                tested=7
            else:
                design=np.column_stack([np.ones(len(time)),time,step,sine,cosine,is_post,tpost]); tested=5
            blocks.append((outcome,design))
        outcome=np.concatenate([block[0] for block in blocks]); design=np.vstack([block[1] for block in blocks])
        try:
            fit=sm.GLM(outcome,design,family=sm.families.Poisson()).fit(scale="X2",maxiter=50,tol=1e-8)
            statistic=fit.params[tested]/fit.bse[tested]; p_value=2*t.sf(abs(statistic),fit.df_resid)
            if not math.isfinite(p_value): raise ValueError("non-finite fit")
            rejected+=p_value<alpha
        except Exception:
            failed+=1
    successful=reps-failed; power=rejected/successful
    return {"power":power,"mcse":math.sqrt(power*(1-power)/successful),"attempted":reps,"successful":successful,"failed":failed,"seed":seed}


its=its_power(inputs["interruptedTimeSeries"]); sw=stepped_wedge_power(inputs["steppedWedge"])
print(f"ITS attainable power {its['power']:.3f} (Monte Carlo SE {its['mcse']:.3f}; {its['successful']}/{its['attempted']} successful fits).")
print(f"Stepped-wedge attainable power {sw['power']:.3f}; size for target {sw['sizeForTarget'] if math.isfinite(sw['sizeForTarget']) else 'not attainable'}.")
print("\n--- HARNESS ---")
for name,value in its.items(): print(f"its_{name}={value:.10f}")
for name,value in sw.items(): print(f"sw_{name}={value:.10f}")
