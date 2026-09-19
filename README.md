# Mallard Code Library

Tested, executable reference implementations of biostatistical analyses in R and Python.

This is the code library behind [Ask Mallard](https://askmallard.com). Ask Mallard turns a
plain-language study description into a full biostatistics plan, and part of that plan is runnable
code. Until now that code was written fresh for every plan and **never executed by anyone**. This
library is the other half: a canonical implementation per method, run against a seeded fixture on
every commit, so the method skeleton a plan builds on is known to work.


## What this library claims, and what it does not

Every entry makes **two separate claims**, and conflating them would defeat the purpose.

**Agreement.** Every executed language produces the same estimates from identical rows, to a
the tolerances recorded in each entry’s `expected.json`. Most outputs use **1e-4**; the GEE working-correlation estimate uses **1e-3** and its model-based standard error uses **5e-3**, with the observed differences and rationale recorded beside those keys.

The first version of this README said 1e-6 and called it "far above optimiser noise". The first
real run refuted that: `survival::clogit` and `statsmodels.ConditionalLogit` agreed to a spread of
5.6e-6 on the log odds ratio and 1.0e-5 on the covariate, so 1e-6 sat *below* the noise floor of
two different optimisers on the same likelihood and failed a correct pair. 1e-4 is about ten times
the observed noise, and an HC0-versus-HC3 variance differs by roughly 4e-3, some forty times the
tolerance, so a real default difference is still caught with room to spare.

That episode is left in this README on purpose. A tolerance chosen by assertion is a check
calibrated against nothing, and the fix was to run the thing and read the number.

This is the claim the library exists for. Ask Mallard's own schema tells the model, on every plan:

> PIN THE DEFAULTS that differ across packages IN THE CODE ITSELF ... silently different defaults
> are how two correct-looking implementations disagree.

That instruction names a real hazard and then asks a language model to be careful. Running the R
and the Python against the same data and comparing the numbers *tests* it. `PROC LOGISTIC` models
the lower ordered value by default, so a 0/1 outcome without `event='1'` inverts every odds ratio
in the output. Reading either file alone would not reveal it. Running both would.

**Recovery.** Those estimates are close to the parameters the fixture actually used. Loose, because
an estimate carries sampling error. It exists to catch implementations that agree with each other
on the same wrong model.

Neither claim implies the other. Agreement alone passes identical mistakes; recovery alone passes a
default mismatch smaller than sampling error. `harness/check.py` reports them separately and refuses
to state the agreement claim at all when fewer than two engines ran.

Installed R and Python package versions are checked against every metadata pin in CI. Execution states in `meta.json` and `expected.json` must agree. Non-finite results, missing numeric outputs and missing recovery targets are failures. Explanatory transformations in `derived_truth` are not separate recovery checks.

### What the first runs actually found

Three real defects in three runs, none of which reading the files would have produced:

1. **A dependency break.** `statsmodels 0.14.2` imports a scipy private helper that newer scipy
   removed. pandas, numpy and statsmodels were pinned; scipy was not, so pip took the latest and
   the import failed outright. Pinning what a file *names* while its dependencies float is not
   pinning.
2. **A CI step asserting its own scope.** The R install list was hardcoded to `survival`, so the
   second entry failed with "there is no package called 'sandwich'" while its `meta.json` had
   declared `sandwich` all along. Package lists are now derived from the entries.
3. **A tolerance that was never measured.** See above.

**What no entry claims is that your analysis is correct.** A working implementation of conditional
logistic regression says nothing about whether conditional logistic regression answers your
question, or whether the right column reached the right argument.

## Verification, stated on every file

| Engine | Status | Why |
|---|---|---|
| R | executed | open source, containerised in CI |
| Python | executed | open source, containerised in CI |

Both languages are executed on the same seeded fixture on every commit, and each file states its
status on its own face. The library ships R and Python only, because those are the two languages Ask
Mallard generates plans in, so those are the two it is built to verify. Earlier revisions also
carried SAS and Stata as unexecuted reference translations; they were removed, because a library
that let a reader assume an unexecuted file was as verified as an executed one would be making
exactly the mistake it exists to prevent.

## Layout

```
lib/<method>/
  meta.json      identity, typed estimand compatibility, Learn page, pinned packages and defaults
  README.md      what it estimates, what it assumes, what it deliberately does not do
  fixture.py     seeded generator, stdlib only
  fixture.csv    COMMITTED, so every language reads byte-identical rows
  expected.json  the two claims and their tolerances
  r.R  python.py
harness/
  check.py       compares the executed engines against each other and against the truth
```

### Why the fixture CSV is committed

The obvious design is to seed each language's own generator. That is quietly wrong: R's Mersenne
Twister stream and NumPy's differ, so "seed 42" produces two different
datasets. Two implementations analysing two datasets cannot be compared at all, and the
cross-language check would be measuring sampling variation while appearing to measure agreement.

The generator is committed too, so the data is reproducible rather than mysterious. Nothing in CI
regenerates it.

### Why expected values are not generated by running this code

They are not recorded from a run of the library's own implementations. A golden regenerated by the
build it is meant to check passes against whatever that build produced, and pins nothing.

Instead the truth is independent of every implementation: the fixture's generating process is
constructed so that the parameters it uses **are** the parameters the model estimates. For
conditional logistic regression, drawing the case within a matched set with probability proportional
to `exp(Xb)` is the conditional likelihood itself.

Where a published worked example with known values exists, that is better still and should be
preferred.

## How Ask Mallard uses it

Plans **compose** from the library rather than reproducing it. The canonical skeleton is supplied to
the code step, and the model adapts variable names, data source and the study's specifics. The
method's pinned defaults survive because nothing rewrote them.

Full mechanical substitution was considered and rejected. Deciding which of a study's columns is the
stratum, the cluster, the period or the competing event type is a judgement, not a find and replace,
and getting it wrong produces code that runs and answers a different question.

## Fixtures are synthetic. Always.

No real dataset, and no patient data of any kind, enters this repository. Every fixture is generated
by seeded code that is committed alongside it.

## Licence

Apache-2.0, with a NOTICE file. This code is meant to be pasted into other people's analyses, so the
licence should not be the thing that stops them; Apache-2.0 stays permissive while adding an explicit
patent grant.
