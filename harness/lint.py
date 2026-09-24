"""Structural lint for what the numeric checks cannot see.

R and Python are both executed on the seeded fixture in CI, so a dropped option usually changes the
numbers and the agreement check fails. `must_appear` is the guard for the cases that DO NOT move a
number -- a pinned tie method, an offset convention, a response level -- and for making sure a
pinned default is set in the CODE rather than merely described in a comment beside it.

WHY THIS IS A SCRIPT AND NOT INLINE CI BASH. The first version was six lines inside the workflow,
and it read a pinned option out of a prose comment -- the unanchored-pattern failure this project
pays for over and over. Inline CI bash is unreviewed code: nothing runs it but the runner, and
nothing tests it at all. As a module it gets controls, below in lint.test.py, including files that
must NOT be flagged.

WHAT THIS LINT ENFORCES. Every expected.json carries a `must_appear` block naming the strings each
language file has to contain. It is enforced against CODE, never comments: a pinned default named in
a comment and absent from the statement below it is exactly the failure the block exists to catch,
and a substring search over the raw file would call that a pass.
"""

import argparse
import json
import re
import sys
from pathlib import Path


def strip_hash_comments(text: str, triple: bool = False) -> str:
    """Remove `#` comments from R or Python while leaving string literals alone.

    A plain `re.sub(r"#.*", "", ...)` would delete half of `sep = "#"` and, worse, would treat a
    `#` inside a quoted string as the start of a comment and strip the rest of the line -- the
    unanchored-pattern failure this repository keeps paying for. So this walks characters and
    tracks quote state instead.

    `triple` is for Python, which R does not have, and it does TWO things: it tracks triple-quoted
    regions so a `#` inside one does not start a comment, and it DROPS their contents. Dropping
    them is the point. A Python docstring is prose in the role a `#` comment plays in R, and the
    first version of this function kept it -- so `must_appear` could declare an option that lived
    only in a module docstring and pass, which is exactly the "described but not set" case the
    check names as the more dangerous of the two. Caught by writing a Python file whose pinned
    option appeared in its docstring and nowhere else.

    A pinned string inside an ordinary single- or double-quoted literal is still CODE and is kept:
    `cov_type='naive'` is an argument, not a description of one.
    """
    out = []
    i, n = 0, len(text)
    quote = None            # the closing delimiter we are waiting for, or None
    drop = False            # inside a triple-quoted region, whose contents are prose
    while i < n:
        ch = text[i]
        if quote:
            if ch == "\\" and i + 1 < n:      # an escaped character cannot close the string
                if not drop:
                    out.append(text[i:i + 2])
                i += 2
                continue
            if text.startswith(quote, i):
                if not drop:
                    out.append(quote)
                i += len(quote)
                quote = None
                drop = False
                continue
            if not drop:
                out.append(ch)
            i += 1
            continue
        if triple and (text[i:i + 3] == '"""' or text[i:i + 3] == "'''"):
            quote = text[i:i + 3]
            drop = True
            i += 3
            continue
        if ch == '"' or ch == "'":
            quote = ch
            out.append(ch)
            i += 1
            continue
        if ch == "#":
            while i < n and text[i] != "\n":
                i += 1
            continue
        out.append(ch)
        i += 1
    return "".join(out)


# filename -> (how to reduce it to code, whether the language cares about case).
# R and Python are both case-sensitive and are matched as written.
LANGUAGES = {
    "r": ("r.R", lambda s: strip_hash_comments(s), True),
    "python": ("python.py", lambda s: strip_hash_comments(s, triple=True), True),
}


def lint_must_appear(entry: Path) -> list:
    """Every string an entry declares in `must_appear` must be present in that file's CODE.

    Two directions, because a rule that only checks what it was handed is a rule asserting its own
    scope -- this repository's most repeated defect, already committed once in the harness's
    hardcoded field names and once in CI's hardcoded R package list:

      * every declared string must appear in the file, and
      * every language file that EXISTS must be declared, or a file could ship with nothing pinned
        and pass in silence.
    """
    expected_path = entry / "expected.json"
    if not expected_path.exists():
        return [f"{entry}: no expected.json"]
    try:
        declared = json.loads(expected_path.read_text()).get("must_appear", {})
    except json.JSONDecodeError as exc:
        return [f"{expected_path}: not valid JSON -- {exc}"]
    if not isinstance(declared, dict):
        return [f"{expected_path}: must_appear is not an object"]

    problems = []
    for lang, (filename, strip, cased) in LANGUAGES.items():
        path = entry / filename
        wanted = declared.get(lang)
        if not path.exists():
            if wanted:
                problems.append(
                    f"{expected_path}: must_appear names {lang} but {filename} does not exist")
            continue
        if not wanted:
            problems.append(
                f"{expected_path}: {filename} exists but must_appear declares nothing for {lang}. "
                f"An implementation with no pinned string is one nothing can check.")
            continue
        raw = path.read_text()
        code = strip(raw)
        haystack, hay_raw = (code, raw) if cased else (code.lower(), raw.lower())
        for needle in wanted:
            probe = needle if cased else needle.lower()
            if probe in haystack:
                continue
            # NAMING WHICH OF THE TWO IT IS, because they are different mistakes: a string that is
            # present only in a comment is a default someone described and did not set, which is
            # the more dangerous of the two and the one a raw substring search would pass.
            where = "present only in a comment" if probe in hay_raw else "absent"
            problems.append(
                f"{path}: must_appear declares {needle!r} and it is {where} in the code")
    return problems


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="lib")
    args = ap.parse_args(argv)

    root = Path(args.root)
    entries = sorted(q.parent for q in root.glob("*/expected.json"))

    problems = []
    for e in entries:
        problems += lint_must_appear(e)

    print(f"entries checked against must_appear: {len(entries)}")
    for p in problems:
        print(f"  FAIL {p}")
    if not problems:
        print("  structural lint passed: every pinned string each entry declared is present in the")
        print("  code rather than only in the prose beside it. Passing here is never a claim that")
        print("  the analysis is right, in either language.")

    # The repository-level checks run from here because CI already calls this file; giving each its
    # own workflow step needs a workflow-scoped change. Each prints its own result.
    import catalogue
    import public_surface
    repo = root.resolve().parent
    stale = catalogue.main(["--readme", str(repo / "README.md"), "--lib", str(root), "--check"])
    surface = public_surface.main(["--repo", str(repo)])
    return 1 if (problems or stale or surface) else 0


if __name__ == "__main__":
    sys.exit(main())
