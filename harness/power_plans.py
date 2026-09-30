"""Validate the analysis-to-power-plan catalogue and expose its resolved records.

Every analysis chooses a typed sizing route. Supported routes point to executable library entries;
unsupported routes still point to executable planning guidance and must explain why an ordinary
power calculation is invalid. This is deliberately a both-directions check: every analysis needs a
route, and every sizing entry must be reached by at least one declared route.
"""

import argparse
import json
from pathlib import Path


ANALYSIS_KIND = "analysis"
SIZING_KIND = "sizing"
FALLBACK_APPROACHES = {
    "precision", "fixed-design", "rule-of-thumb", "bespoke-simulation", "not-applicable"
}
ROUTE_APPROACHES = FALLBACK_APPROACHES | {"analytic-power"}


def read_json(path: Path):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"{path}: could not read valid JSON -- {exc}") from exc


def resolved_plans(repo: Path):
    routes_doc = read_json(repo / "sizing-routes.json")
    routes = routes_doc.get("routes")
    fallbacks = routes_doc.get("fallbacks")
    if not isinstance(routes, dict) or not routes or not isinstance(fallbacks, dict):
        raise ValueError("sizing-routes.json: routes and fallbacks must be non-empty objects")
    if set(fallbacks) != FALLBACK_APPROACHES:
        raise ValueError("sizing-routes.json: fallback approaches do not match the public vocabulary")

    metas = {}
    for path in sorted((repo / "lib").glob("*/meta.json")):
        meta = read_json(path)
        if meta.get("id") != path.parent.name:
            raise ValueError(f"{path}: id does not equal its directory name")
        metas[meta["id"]] = (path, meta)

    analyses = {key: value for key, value in metas.items()
                if value[1].get("kind", ANALYSIS_KIND) == ANALYSIS_KIND}
    sizing_entries = {key: value for key, value in metas.items()
                      if value[1].get("kind", ANALYSIS_KIND) == SIZING_KIND}
    unknown = sorted((path, meta.get("kind")) for path, meta in metas.values()
                     if meta.get("kind", ANALYSIS_KIND) not in {ANALYSIS_KIND, SIZING_KIND})
    if unknown:
        raise ValueError(f"{unknown[0][0]}: unknown entry kind {unknown[0][1]!r}")
    if not analyses:
        raise ValueError("lib/: no analysis entries found")

    def validate_route(route_key, route, required_inputs=True):
        if not isinstance(route, dict):
            raise ValueError(f"sizing-routes.json: {route_key} is not an object")
        approach = route.get("approach", route_key if route_key in FALLBACK_APPROACHES else None)
        if approach not in ROUTE_APPROACHES:
            raise ValueError(f"sizing-routes.json: {route_key} has invalid approach {approach!r}")
        power_id = route.get("entry")
        target = sizing_entries.get(power_id)
        if not target:
            raise ValueError(f"sizing-routes.json: {route_key} names missing sizing entry {power_id!r}")
        supported = target[1].get("sizing_calcs")
        if not isinstance(supported, list) or route_key not in supported:
            raise ValueError(f"{target[0]}: sizing_calcs does not cover declared route {route_key!r}")
        if not isinstance(route.get("learn"), str) or not route["learn"]:
            raise ValueError(f"sizing-routes.json: {route_key} has no Learn page")
        required = route.get("requiredInputs")
        if required_inputs and (not isinstance(required, list) or not required
                                or any(not isinstance(v, str) or not v for v in required)):
            raise ValueError(f"sizing-routes.json: {route_key} has no usable requiredInputs")
        if not isinstance(route.get("nextStep"), str) or len(route["nextStep"].strip()) < 30:
            raise ValueError(f"sizing-routes.json: {route_key} has no substantive nextStep")
        return route

    # Validate the complete public vocabulary, not only the routes the current 86 analyses happen
    # to use. Otherwise an unused future route can silently point to no code or no Learn page.
    for calc, route in routes.items():
        validate_route(calc, route)
    for approach, route in fallbacks.items():
        validate_route(approach, route)
    declared_targets = {route["entry"] for route in routes.values()} | {
        route["entry"] for route in fallbacks.values()
    }
    orphaned = sorted(set(sizing_entries) - declared_targets)
    if orphaned:
        raise ValueError(f"sizing entries reached by no declared route: {', '.join(orphaned)}")

    out = {}
    for entry_id, (path, meta) in analyses.items():
        sizing = meta.get("sizing")
        if not isinstance(sizing, dict):
            raise ValueError(f"{path}: every analysis must declare a sizing object")
        calc = sizing.get("calc")
        approach = sizing.get("approach")
        if not isinstance(calc, str) or not calc:
            raise ValueError(f"{path}: sizing.calc must be a non-empty string")

        if calc == "none":
            route = fallbacks.get(approach)
            why = sizing.get("why")
            if not route:
                raise ValueError(f"{path}: calc none needs a typed fallback approach")
            if not isinstance(why, str) or len(why.strip()) < 40:
                raise ValueError(f"{path}: calc none needs a substantive sizing.why")
        else:
            route = routes.get(calc)
            why = sizing.get("why") or sizing.get("note") or ""
            if not route:
                raise ValueError(f"{path}: sizing.calc {calc!r} has no canonical route")
            if approach != route.get("approach"):
                raise ValueError(
                    f"{path}: sizing.approach {approach!r} disagrees with route {route.get('approach')!r}")

        power_id = route["entry"]
        required = route["requiredInputs"]
        out[entry_id] = {
            "calc": calc,
            "approach": approach,
            "why": why,
            "requiredInputs": required,
            "nextStep": route["nextStep"],
            "learn": route["learn"],
            "entry": power_id,
        }

    return out


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", default=str(Path(__file__).resolve().parent.parent))
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args(argv)
    plans = resolved_plans(Path(args.repo))
    if args.json:
        print(json.dumps(plans, indent=2, sort_keys=True))
    else:
        counts = {}
        for value in plans.values():
            counts[value["approach"]] = counts.get(value["approach"], 0) + 1
        print(f"analysis sizing routes checked: {len(plans)}")
        print("  " + ", ".join(f"{key}={counts[key]}" for key in sorted(counts)))
        print("  every analysis has code, Learn guidance and a generated-plan route")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
