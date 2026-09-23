# An incidence rate (events per person-year) with an exact Poisson 95% interval
#
# The Python equivalent of r.R. statsmodels' confint_poisson is an implementation independent of R's
# poisson.test, so the two agreeing is evidence rather than one formula typed twice.

import pandas as pd
from statsmodels.stats.rates import confint_poisson

d = pd.read_csv("fixture.csv")
assert not d.isna().any().any()
assert (d.followup_years > 0).all() and (d.events >= 0).all()

events = int(d.events.sum())
person_years = float(d.followup_years.sum())
rate = events / person_years

# PINNED: method="exact-c", the central exact (Garwood) interval, which is what R's poisson.test
# reports. confint_poisson offers about a dozen methods, and "exact-c" is not the only one whose
# name sounds exact. The exposure argument is the person-time; leaving it at 1 gives the interval
# for the count.
rate_lcl, rate_ucl = confint_poisson(events, person_years, method="exact-c", alpha=0.05)

print(f"{events} events over {person_years:.1f} person-years")
print(f"rate {rate:.4f} per person-year ({1000 * rate:.1f} per 1000), "
      f"exact 95% {rate_lcl:.4f} to {rate_ucl:.4f}")

print("\n--- HARNESS ---")
print(f"rate={rate:.10f}")
print(f"rate_lcl={float(rate_lcl):.10f}")
print(f"rate_ucl={float(rate_ucl):.10f}")
print(f"events={events}")
print(f"person_years={person_years:.10f}")
print(f"n={len(d)}")
