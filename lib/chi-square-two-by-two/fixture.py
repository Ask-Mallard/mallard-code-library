"""Seeded fixture for a two-by-two table: chi-square, Fisher's exact test and the odds ratio.

The outcome is drawn from a Bernoulli distribution whose probability depends only on the exposure,
so the odds ratio the generator implies IS the parameter the table estimates -- no model, no link
function, nothing to misspecify.

EXPECTED COUNTS ARE COMFORTABLY LARGE HERE, AND THAT IS DELIBERATE. Yates' continuity correction is
defended as a small-sample device, so a fixture with sparse cells would let a reader dismiss the
disagreement between packages as an edge case. It is not one: with two thousand patients and every
expected count in the hundreds, the corrected and uncorrected chi-square still differ by 0.715,
because R and Python apply Yates' correction by default and every file here switches it off.

Stdlib only.
"""

import csv
import random
from pathlib import Path

SEED = 20260908
N = 2000
P_EXPOSED = 0.5
RISK_EXPOSED = 0.35
RISK_UNEXPOSED = 0.20
# true odds ratio (0.35/0.65) / (0.20/0.80) = 2.153846..., log 0.767255...


def main() -> None:
    rng = random.Random(SEED)
    rows = []
    for i in range(1, N + 1):
        exposed = 1 if rng.random() < P_EXPOSED else 0
        risk = RISK_EXPOSED if exposed else RISK_UNEXPOSED
        outcome = 1 if rng.random() < risk else 0
        rows.append({"id": i, "exposed": exposed, "outcome": outcome})

    out = Path(__file__).with_name("fixture.csv")
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["id", "exposed", "outcome"])
        w.writeheader()
        w.writerows(rows)

    a = sum(1 for r in rows if r["exposed"] == 1 and r["outcome"] == 1)
    b = sum(1 for r in rows if r["exposed"] == 1 and r["outcome"] == 0)
    c = sum(1 for r in rows if r["exposed"] == 0 and r["outcome"] == 1)
    d = sum(1 for r in rows if r["exposed"] == 0 and r["outcome"] == 0)
    print(f"wrote {out.name}: {len(rows)} rows")
    print(f"            outcome=1  outcome=0")
    print(f"  exposed    {a:8d}   {b:8d}")
    print(f"  unexposed  {c:8d}   {d:8d}")
    n = a + b + c + d
    for lab, r_, c_ in (("a", a, (a + b) * (a + c)), ("b", b, (a + b) * (b + d)),
                        ("c", c, (c + d) * (a + c)), ("d", d, (c + d) * (b + d))):
        print(f"  expected count for {lab}: {c_ / n:.1f}")
    print(f"  crude risks: exposed {a / (a + b):.4f} (true {RISK_EXPOSED}), "
          f"unexposed {c / (c + d):.4f} (true {RISK_UNEXPOSED})")
    print(f"  sample odds ratio {(a * d) / (b * c):.6f}   true "
          f"{(RISK_EXPOSED / (1 - RISK_EXPOSED)) / (RISK_UNEXPOSED / (1 - RISK_UNEXPOSED)):.6f}")


if __name__ == "__main__":
    main()
