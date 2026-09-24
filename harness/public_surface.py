"""Keep development-process text out of a public repository.

This library is public and is read by people who have never seen the project that maintains it. The
text a maintainer's tooling produces while building it (decision labels, internal file names,
session links, batch codes) means nothing to them and says more about the maintainer's process than
the library should. Every pattern below was found in this repository's tracked files or commit
messages before the check existed.

The patterns are ANCHORED to the phrase, never a bare word: `spec` is also a variable holding a
module spec and a specificity, `owner` is in the Apache licence, and `batch` is ordinary English.
lint.test.py pins each of those as a line that must NOT be flagged.

Usage:
    python harness/public_surface.py                 # scan every tracked file
    python harness/public_surface.py --text body.md  # scan one file, e.g. a pull request body
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

# (pattern, why it does not belong in a public file)
FORBIDDEN = [
    (r"\bowner decisions?\b", "an internal decision label; state the rule, not who decided it"),
    (r"\bCLAUDE\.md\b", "an internal instructions file; cite the published source for the rule"),
    (r"\bdocs/specs/", "a path in a private repository"),
    (r"\bschema prompt\b", "describes a private product's prompt"),
    (r"\bevaluation bank\b", "a private evaluation set"),
    (r"\bClaude-Session:", "a link to a private working session"),
    (r"claude\.ai/code/session_", "a link to a private working session"),
    (r"\bbatch [A-I]-\d\b", "an internal build-batch code; describe the entries instead"),
    (r"\bcategory [A-I]\b", "an internal build-plan code; describe the entries instead"),
    (r"\bStatistical correctness\b", "a heading in a private instructions file"),
    (r"\bAgent's transcription\b", "a working-session transcript"),
]
COMPILED = [(re.compile(p), why) for p, why in FORBIDDEN]

# This file and its controls necessarily spell out the patterns they look for.
EXEMPT = {"harness/public_surface.py", "harness/lint.test.py"}


def scan_text(text, label="<text>"):
    problems = []
    for lineno, line in enumerate(text.splitlines(), 1):
        for rx, why in COMPILED:
            m = rx.search(line)
            if m:
                problems.append(f"{label}:{lineno}: {m.group(0)!r} -- {why}")
    return problems


def tracked_files(repo):
    """Every tracked file, from git. A directory walk would also scan ignored local output nobody
    sees on GitHub, so a missing git is an error rather than a fallback."""
    out = subprocess.run(["git", "ls-files"], cwd=repo, capture_output=True, text=True, check=True)
    return [f for f in out.stdout.splitlines() if f]


def scan_repo(repo):
    repo = Path(repo)
    problems, scanned = [], 0
    for rel in tracked_files(repo):
        if rel in EXEMPT:
            continue
        try:
            text = (repo / rel).read_text(encoding="utf-8")
        except (UnicodeDecodeError, FileNotFoundError):
            continue
        scanned += 1
        problems += scan_text(text, rel)
    if scanned == 0:
        problems.append("no tracked files were scanned; the check would pass on anything")
    return problems, scanned


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default=str(Path(__file__).resolve().parent.parent))
    ap.add_argument("--text", help="scan this one file instead of the repository")
    args = ap.parse_args(argv)
    if args.text:
        problems = scan_text(Path(args.text).read_text(encoding="utf-8"), args.text)
        print(f"public-surface check of {args.text}")
    else:
        problems, n = scan_repo(args.repo)
        print(f"public-surface check: {n} tracked files scanned")
    for p in problems:
        print(f"  FAIL {p}")
    if not problems:
        print("  no development-process text found")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
