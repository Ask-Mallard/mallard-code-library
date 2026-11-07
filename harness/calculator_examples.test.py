"""Run calculator entries against each other and saved R reference outputs.

The simulated-parameter recovery checks remain in force. Saved R results add a
fixed-fixture regression check, including statistics and intervals that have no
population parameter in an entry's recovery block. References are recorded only
from executed R output; this test never generates or updates expected values.
"""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

from check import compare, parse_harness_block

ROOT = Path(__file__).resolve().parent.parent


def validate_reference(entry, expected, meta):
    reference = expected.get("r_reference")
    if not reference:
        raise ValueError(f"{entry.name}: R reference pending")
    if reference["r_version"] != "4.6.1":
        raise ValueError(f"{entry.name}: reference requires R 4.6.1")
    if reference["snapshot"] != "https://packagemanager.posit.co/cran/__linux__/noble/2026-09-01":
        raise ValueError(f"{entry.name}: unexpected package snapshot")
    if reference["packages"] != meta["packages"]["r"]:
        raise ValueError(f"{entry.name}: reference package versions differ from metadata")
    digest = hashlib.sha256((entry / expected["fixture"]).read_bytes()).hexdigest()
    if reference["fixture_sha256"] != digest:
        raise ValueError(f"{entry.name}: reference fixture digest differs")
    values = reference["values"]
    if not values or any(type(v) not in (int, float) for v in values.values()):
        raise ValueError(f"{entry.name}: reference must contain numeric results")
    return values


def main():
    coverage = json.loads((ROOT / "harness/calculator_coverage.json").read_text())
    entries = sorted({name for group in coverage.values() for name in group})
    assert len(coverage) == 17, "calculator category missing from coverage inventory"
    assert all(coverage.values()), "calculator category without a library entry"
    failures = []
    with tempfile.TemporaryDirectory() as directory:
        logs = Path(directory)
        for name in entries:
            entry = ROOT / "lib" / name
            try:
                expected = json.loads((entry / "expected.json").read_text())
                meta = json.loads((entry / "meta.json").read_text())
                values = validate_reference(entry, expected, meta)
                version = subprocess.check_output(
                    ["Rscript", "-e", "cat(as.character(getRversion()))"], text=True).strip()
                assert version == "4.6.1", f"R version {version} differs from reference"
                # Package installation can report success with missing packages. Check actual
                # versions before using an engine, including R's bundled recommended packages.
                for package, wanted in meta["packages"]["r"].items():
                    actual = subprocess.check_output(
                        ["Rscript", "-e", f'cat(as.character(packageVersion("{package}")))'],
                        text=True).strip()
                    assert actual == wanted, f"{package}: installed {actual}, required {wanted}"
                results = {}
                arguments = [sys.executable, str(ROOT / "harness/check.py"), "--entry", str(entry)]
                for language, command in (
                    ("r", ["Rscript", "r.R"]),
                    ("python", [sys.executable, "python.py"]),
                ):
                    output = subprocess.check_output(command, cwd=entry, text=True)
                    results[language] = parse_harness_block(output)
                    assert set(results[language]) == set(values), (
                        f"{language}: reference and output keys differ")
                    path = logs / f"{name}.{language}.txt"
                    path.write_text(output)
                    arguments.extend(["--output", f"{language}={path}"])
                # Execute the normal recovery/agreement check, then use the unchanged
                # comparator with R's recorded values as a third set of results.
                subprocess.run(arguments, check=True)
                findings = compare(
                    {**results, "r_reference": values},
                    {"agreement": expected["agreement"], "truth": {}},
                )
                assert not findings, "\n".join(findings)
                print(f"R_REFERENCE|{name}|r,python|TESTED|pass", flush=True)
            except (AssertionError, KeyError, ValueError, OSError, subprocess.SubprocessError) as exc:
                failures.append(f"{name}: {exc}")
    for failure in failures:
        print(f"FAIL {failure}", flush=True)
    print(f"calculator entries: {len(entries)}; passed: {len(entries) - len(failures)}; "
          f"failed or pending: {len(failures)}", flush=True)
    return bool(failures)


if __name__ == "__main__":
    sys.exit(main())
