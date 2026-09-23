# Statistical process control: u-, p- and individuals charts

Control charts show whether a process has changed. Pick the chart for the data:

| Data | Chart | This fixture |
|---|---|---|
| counts over varying exposure | **u-chart** | central-line infections per 1,000 catheter-days: rises from 2 to 5 at month 19; signals from month 19 |
| proportions over varying denominators | **p-chart** | 30-day readmissions: steady at 0.12; no signal |
| one continuous value per period | **individuals (XmR)** | mean length of stay: steady at 5.0; no signal |

## The defaults this entry pins

- **Limits from a baseline, then applied forward**: months 1 to 12 set the centre and 3-sigma limits for
  all 24 months. Limits recomputed from data that include a change absorb it.
- u- and p-chart limits vary with each month's denominator; the individuals chart's sigma is the average
  moving range divided by 1.128.
- A signal is a point beyond a limit (run rules not applied).

## What the simulations show

- **Rare events hide changes.** At 2 to 3 infections a month, a doubling of the rate did not signal; the
  fixture uses a larger hospital.
- **Short baselines give noisy limits**: 12 baseline points produced about 1% false alarms per point,
  against the nominal 0.27%. Use 20 to 25 points where possible.
- A signal shows a process changed; attributing the change to an intervention needs a design (see
  `interrupted-time-series-segmented`).

## Verification

| Engine | Status |
|---|---|
| R (qcc 2.7) | executed in CI |
| Python (written out) | executed in CI |
