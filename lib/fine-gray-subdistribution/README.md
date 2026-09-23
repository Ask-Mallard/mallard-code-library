# Fine-Gray regression for the cumulative incidence of one cause

When a competing event (death) can prevent the event of interest (relapse), the question "does
treatment change how many people relapse?" is about the **cumulative incidence** of relapse. The
Fine-Gray model answers it with a **subdistribution hazard ratio**.

## How to read it

The subdistribution risk set keeps people who already died of the competing cause, so the ratio is not
a rate among people at risk. Read it through the cumulative incidence it implies ("treatment lowers the
proportion who relapse by 3 years"), and report the cumulative incidence curves beside it
(`competing-risks-cumulative-incidence`). For the RATE of relapse among those still alive, use
`cause-specific-cox`; many analyses report both.

| Model | Log hazard ratio here | Question |
|---|---|---|
| **Fine-Gray** | −0.625 (truth −0.5) | effect on the cumulative incidence of relapse |
| cause-specific Cox | −0.532 | effect on the relapse rate among those event-free |

## The defaults this entry pins

`failcode = 1` and `cencode = 0` are named even though they are `crr`'s defaults: data coding death as
1 would silently model death. `cov1` must be a numeric matrix; `crr` takes no formula.

## R only

No standard Python package implements Fine-Gray regression, and hand-writing the censoring-weighted
risk-set expansion is exactly what should come from a maintained implementation. The recovery claim is
checked by fitting 100 generated fixtures in R, and a second R implementation (`survival::finegray`
with a weighted, clustered `coxph`) agrees with `crr` to 1e-5.

## The fixture

1000 people, generated with Fine and Gray's own 1999 simulation design, so the subdistribution hazard
ratio is known exactly (log −0.5). 262 relapses, 486 competing deaths.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | not applicable (see above) |
