# Two-period crossover trial (AB/BA)

Each patient receives both treatments in random order (AB or BA) with a washout between. Each patient's
period-1 minus period-2 difference contains the treatment effect, with opposite signs in the two
sequences, and the period effect, with the same sign. Comparing those differences between sequences
cancels the period effect:

treatment effect (B − A) = (mean difference in BA − mean difference in AB) / 2

## The trap

A paired comparison of B with A that ignores the period is biased whenever the sequences are unequal in
size: the period effect leaks in by period effect × (n_AB − n_BA) / n. With a period effect of 4 and
35 against 25 patients, that is 0.67. On this fixture:

| Analysis | Estimate (truth −3) |
|---|---|
| **period-adjusted** | −2.55 |
| paired B − A ignoring period | −1.98 |

## Carryover

The analysis assumes no carryover, which the washout must make plausible. A test for carryover has little
power and should not decide the analysis.

## The fixture

60 patients (35 AB, 25 BA); patient level SD 6; period 2 adds 4; B − A = −3; noise SD 3.

## Verification

| Engine | Status |
|---|---|
| R | executed in CI |
| Python | executed in CI |
