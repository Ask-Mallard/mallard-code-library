"""Reproduce the committed long-form fixed-design sizing scenarios with the standard library only."""

import csv
import io
from pathlib import Path


FIXTURE = """calc,key,value
interruptedTimeSeries,preTimes,38
interruptedTimeSeries,postTimes,36
interruptedTimeSeries,series,2
interruptedTimeSeries,baselineEvents,6.75
interruptedTimeSeries,controlBaselineEvents,4.9
interruptedTimeSeries,rateRatio,0.7
interruptedTimeSeries,slopeRatio,1
interruptedTimeSeries,nbSize,30
interruptedTimeSeries,stepAt,25
interruptedTimeSeries,stepRatio,0.6
interruptedTimeSeries,seasonalAmplitude,0.15
interruptedTimeSeries,alpha,0.05
interruptedTimeSeries,reps,500
interruptedTimeSeries,seed,20260915
steppedWedge,clusters,8
steppedWedge,periods,5
steppedWedge,clusterPeriodSize,30
steppedWedge,icc,0.03
steppedWedge,cac,1
steppedWedge,baselineRate,0.2
steppedWedge,targetRate,0.15
steppedWedge,alpha,0.05
steppedWedge,power,0.8
"""


def rows():
    yield from csv.DictReader(io.StringIO(FIXTURE))


if __name__ == "__main__":
    Path(__file__).with_name("fixture.csv").write_text(FIXTURE, encoding="utf-8")
