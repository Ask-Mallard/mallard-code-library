# Standard power and precision planning

`r.R` and `python.py` implement the same typed dispatcher used by a generated Mallard power plan.
The analysis entry chooses the `calc`; the code never chooses a method from prose. The fixture runs
every supported route with synthetic assumptions.

The returned count is the analyzable count. Attrition, missingness and eligibility losses are
separate visible assumptions. `logisticEPV` and `linearSPV` return floors and are not power
calculations. Precision routes return the sample needed for the requested interval width. Unsupported
or underspecified methods use the separate `power-plan-guidance` entry instead of borrowing one of
these formulas.
