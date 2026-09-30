"""Failure controls for the catalogue-wide power-plan invariant."""

import json
from pathlib import Path
import shutil
import tempfile

from power_plans import resolved_plans


SOURCE = Path(__file__).resolve().parent.parent


def copied_repo():
    tmp = tempfile.TemporaryDirectory()
    target = Path(tmp.name) / "repo"
    shutil.copytree(SOURCE, target, ignore=shutil.ignore_patterns(".git", "__pycache__"))
    return tmp, target


def rejected(label, mutate):
    tmp, repo = copied_repo()
    try:
        mutate(repo)
        try:
            resolved_plans(repo)
        except ValueError:
            print(f"PASS {label} rejected")
        else:
            raise AssertionError(f"{label} falsely passed")
    finally:
        tmp.cleanup()


plans = resolved_plans(SOURCE)
analysis_count = sum(1 for path in (SOURCE / "lib").glob("*/meta.json")
                     if json.loads(path.read_text()).get("kind", "analysis") == "analysis")
assert len(plans) == analysis_count == 86
routes_doc = json.loads((SOURCE / "sizing-routes.json").read_text())
declared_approaches = {route["approach"] for route in routes_doc["routes"].values()} | set(routes_doc["fallbacks"])
assert declared_approaches == {
    "analytic-power", "precision", "fixed-design", "rule-of-thumb", "bespoke-simulation",
    "not-applicable",
}


def drop_sizing(repo):
    path = repo / "lib" / "welch-two-sample-t-test" / "meta.json"
    value = json.loads(path.read_text())
    del value["sizing"]
    path.write_text(json.dumps(value))


def add_future_without_sizing(repo):
    source = repo / "lib" / "welch-two-sample-t-test"
    target = repo / "lib" / "future-analysis"
    shutil.copytree(source, target)
    meta_path = target / "meta.json"
    value = json.loads(meta_path.read_text())
    value["id"] = "future-analysis"
    del value["sizing"]
    meta_path.write_text(json.dumps(value))


def strip_fallback_reason(repo):
    path = repo / "lib" / "bland-altman-limits" / "meta.json"
    value = json.loads(path.read_text())
    value["sizing"]["why"] = ""
    path.write_text(json.dumps(value))


def lie_about_approach(repo):
    path = repo / "lib" / "welch-two-sample-t-test" / "meta.json"
    value = json.loads(path.read_text())
    value["sizing"]["approach"] = "rule-of-thumb"
    path.write_text(json.dumps(value))


def remove_power_code(repo):
    shutil.rmtree(repo / "lib" / "power-plan-core")


def drop_learn(repo):
    path = repo / "sizing-routes.json"
    value = json.loads(path.read_text())
    del value["routes"]["twoMeans"]["learn"]
    path.write_text(json.dumps(value))


def drop_fallback_inputs(repo):
    path = repo / "sizing-routes.json"
    value = json.loads(path.read_text())
    del value["fallbacks"]["bespoke-simulation"]["requiredInputs"]
    path.write_text(json.dumps(value))


def leave_future_route_without_code(repo):
    path = repo / "sizing-routes.json"
    value = json.loads(path.read_text())
    value["routes"]["futureCalc"] = {
        "approach": "analytic-power",
        "entry": "power-plan-core",
        "learn": "concept.power",
        "requiredInputs": ["effect"],
        "nextStep": "Implement and validate the future calculation before recommending a sample size."
    }
    path.write_text(json.dumps(value))


for label, mutation in [
    ("existing analysis without sizing", drop_sizing),
    ("future analysis without sizing", add_future_without_sizing),
    ("fallback without explanation", strip_fallback_reason),
    ("analysis contradicting its route", lie_about_approach),
    ("route without public code", remove_power_code),
    ("route without Learn", drop_learn),
    ("fallback without required inputs", drop_fallback_inputs),
    ("future route absent from sizing code", leave_future_route_without_code),
]:
    rejected(label, mutation)

print("PASS all 86 analyses and future-entry controls")
