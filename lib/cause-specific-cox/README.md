# Cause-specific Cox regression under competing risks

When another event (death) can prevent the event of interest (relapse), the cause-specific hazard
ratio compares the **rate** of relapse among people still alive and relapse-free. Deaths are censored
**for this model only**.

## Rate, not risk

The cause-specific hazard ratio does not say how many people relapse by 5 years: that also depends on
how fast death removes people. Report the cumulative incidence (`competing-risks-cumulative-incidence`)
beside it, and use a Fine-Gray model when the question is the effect on that cumulative incidence. A
Kaplan-Meier curve that censors deaths overstates the cumulative incidence of relapse.

## The default this entry pins

| Approach | Log hazard ratio here | What it estimates |
|---|---|---|
| **event = relapse only** | −0.378 | the cause-specific hazard ratio (truth −0.357) |
| relapse or death as the event | −0.243 | a composite, a different question |
| dropping people who died | (biased) | removes person-time that was at risk of relapse |

## The fixture

2500 people; relapse hazard 0.10 a year (× 0.7 if treated), death 0.06 a year, follow-up up to 5 years.
658 relapses, 484 deaths.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
