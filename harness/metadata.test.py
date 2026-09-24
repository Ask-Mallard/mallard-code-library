"""Failure controls against real entry metadata, in disposable copies."""
import json
from pathlib import Path
import shutil
import tempfile
from metadata import declarations

root = Path(__file__).resolve().parent.parent / "lib"
declarations(root)  # Byte-faithful negative control: the whole committed library.
cases = ("missing-state", "contradictory-state", "conflicting-version", "missing-source", "missing-compatibility", "unknown-target", "unknown-framing")
for case in cases:
    with tempfile.TemporaryDirectory() as tmp:
        target = Path(tmp) / "lib"
        shutil.copytree(root, target)
        path = target / "cox-proportional-hazards/meta.json"
        value = json.loads(path.read_text())
        if case == "missing-state":
            del value["engines"]["r"]
        elif case == "contradictory-state":
            value["engines"]["r"] = "not-executed"
        elif case == "conflicting-version":
            value["packages"]["r"]["survival"] = "0.0"
        elif case == "missing-source":
            path.with_name("r.R").unlink()
        elif case == "unknown-target":
            value["compatibility"]["target"] = ["prediction"]
        elif case == "unknown-framing":
            value["compatibility"]["framing"] = ["noninferior"]
        else:
            del value["compatibility"]["effectScale"]
        path.write_text(json.dumps(value))
        try:
            declarations(target)
        except ValueError:
            print(f"PASS {case} rejected")
        else:
            raise AssertionError(f"{case} falsely passed")
print("PASS actual library declarations")
