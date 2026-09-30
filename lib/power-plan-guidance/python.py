"""Structured planning guidance when an ordinary power calculation is invalid."""

import csv


def planning_guidance(approach, why, required_inputs, next_step):
    allowed = {"precision", "fixed-design", "rule-of-thumb", "bespoke-simulation", "not-applicable"}
    if approach not in allowed or not why or not required_inputs or not next_step:
        raise ValueError("a typed approach, reason, required inputs and next step are required")
    return {
        "approach": approach,
        "powerCalculated": False,
        "sampleSizeCalculated": False,
        "why": why,
        "requiredInputs": required_inputs.split(";"),
        "nextStep": next_step,
    }


with open("fixture.csv", newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle))
plans = [planning_guidance(row["approach"], row["why"], row["required_inputs"], row["next_step"]) for row in rows]
for plan in plans:
    print(f"{plan['approach']}\n  Why: {plan['why']}\n  Required: {', '.join(plan['requiredInputs'])}\n  Next: {plan['nextStep']}")
print("\n--- HARNESS ---")
print(f"route_count={len(plans)}")
print(f"numeric_sample_sizes={sum(plan['sampleSizeCalculated'] for plan in plans)}")
print(f"routes_with_next_step={sum(bool(plan['nextStep']) for plan in plans)}")
