"""Validate declared execution states and enforce the installed package versions.

Version strings are claims about an executed environment, not decoration. The dated
R repository previously installed versions different from every declared pin but one.
"""
import argparse
import importlib.metadata
import json
from pathlib import Path
import re
import subprocess

FILES = {"r": "r.R", "python": "python.py"}
COMPATIBILITY_KEYS = {"target", "outcomeType", "summaryMeasure", "effectScale", "samplingStructure", "competingEvents", "repeatedMeasures", "missingData", "variant"}


def declarations(root):
    packages = {"r": {}, "python": {}}
    for path in sorted(Path(root).glob("*/meta.json")):
        meta = json.loads(path.read_text())
        expected = json.loads(path.with_name("expected.json").read_text())
        compatibility = meta.get("compatibility")
        if not isinstance(compatibility, dict) or set(compatibility) != COMPATIBILITY_KEYS:
            raise ValueError(f"{path}: incomplete typed compatibility profile")
        for key, values in compatibility.items():
            if not isinstance(values, list) or not values or any(type(v) is not (bool if key in ("competingEvents", "repeatedMeasures") else str) for v in values):
                raise ValueError(f"{path}: invalid compatibility values for {key}")
        states = meta.get("engines", {})
        for lang, filename in FILES.items():
            state = states.get(lang)
            if state not in ("executed", "not-executed", "not-applicable"):
                raise ValueError(f"{path}: missing or invalid {lang} execution state")
            if path.with_name(filename).exists() != (state != "not-applicable"):
                raise ValueError(f"{path}: {lang} file contradicts execution state")
            if expected.get("engines", {}).get(lang, {}).get("executed") is not (state == "executed"):
                raise ValueError(f"{path}: {lang} execution declarations disagree")
        for lang in packages:
            for name, version in meta.get("packages", {}).get(lang, {}).items():
                if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_.-]*", name) or not isinstance(version, str):
                    raise ValueError(f"{path}: invalid package declaration")
                previous = packages[lang].get(name)
                if previous is not None and previous != version:
                    raise ValueError(f"{path}: conflicting {lang} pin for {name}: {previous} / {version}")
                packages[lang][name] = version
    if not any(packages.values()):
        raise ValueError("no package declarations found")
    return packages


def check_installed(packages):
    actual = {"r": {}, "python": {}}
    for name in packages["r"]:
        value = subprocess.check_output(["Rscript", "-e", f'cat(packageDescription("{name}")$Version)'], text=True).strip()
        actual["r"][name] = value
    for name in packages["python"]:
        actual["python"][name] = importlib.metadata.version(name)
    for lang, pins in packages.items():
        for name, version in pins.items():
            if actual[lang][name] != version:
                raise ValueError(f"{lang} {name}: installed {actual[lang][name]}, declared {version}")
    return actual


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="lib")
    ap.add_argument("--installed", action="store_true")
    args = ap.parse_args()
    pins = declarations(args.root)
    print(json.dumps(check_installed(pins) if args.installed else pins, indent=2))
