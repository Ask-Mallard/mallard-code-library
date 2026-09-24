"""Controls for lint.py. Stdlib only; run with `python harness/lint.test.py`.

Every rule gets BOTH directions: a file it must flag, and a file it must not. The false positive
this lint produced on its first CI run -- a `#` read out of a string literal -- is pinned here as a
negative control, built from the exact text that tripped it rather than a paraphrase. A paraphrased
fixture has passed in this project's history while proving nothing.
"""

import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import json  # noqa: E402
from lint import lint_must_appear, strip_hash_comments  # noqa: E402

FAILURES = []


def check(name, condition, detail=""):
    if condition:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name} {detail}")
        FAILURES.append(name)


def write(tmp, name, text):
    p = Path(tmp) / name
    p.write_text(text)
    return p


def main():
    with tempfile.TemporaryDirectory() as tmp:
        print("comment stripping (R and Python)")

        # A `#` inside a string is not a comment. A regex-based stripper deletes the rest of the
        # line here, and `must_appear` would then report a string it can plainly see as absent.
        r_hash = 'sep <- "#"\nfit <- coxph(Surv(time, event) ~ exposed, ties = "efron")\n'
        check("a # inside an R string is not treated as a comment",
              'ties = "efron"' in strip_hash_comments(r_hash))
        check("a real R comment IS stripped",
              "efron" not in strip_hash_comments('fit <- lm(y ~ x)  # ties = "efron"\n'))

        # A PYTHON DOCSTRING IS PROSE, and is dropped. This test asserted the opposite until
        # 2026-09-08, when a python.py declared cov_struct=Exchangeable() in must_appear, wrote it
        # in its module docstring, called the constructor a different way in the code, and PASSED.
        # That is the "described but not set" case the check names as the more dangerous of the
        # two -- available in the one language whose convention is to describe things in a string.
        py_doc = '''"""Docstring mentioning # and cov_type="HC3" in prose."""
res = mod.fit(cov_type="HC1")
'''
        stripped = strip_hash_comments(py_doc, triple=True)
        check("a Python docstring's contents are DROPPED, like the comment they are",
              'cov_type="HC3"' not in stripped)
        check("code after a docstring containing a # survives it",
              'cov_type="HC1"' in stripped)
        check("an ordinary Python string literal is CODE and is kept",
              "cov_type='naive'" in strip_hash_comments(
                  "se = fit.standard_errors(cov_type='naive')\n", triple=True))
        check("a real Python comment IS stripped",
              "HC3" not in strip_hash_comments('res = mod.fit()  # HC3 would go here\n',
                                               triple=True))
        # BOTH DIRECTIONS on the real case: the same file, the option in the docstring only,
        # against the same file with it in the call.
        only_prose = '''"""Uses cov_struct=Exchangeable() throughout."""
model = sm.GEE(y, X, groups=g, cov_struct=Independence())
'''
        in_code = '''"""Uses an exchangeable working correlation."""
model = sm.GEE(y, X, groups=g, cov_struct=Exchangeable())
'''
        check("an option present only in the docstring is NOT code",
              "cov_struct=Exchangeable()" not in strip_hash_comments(only_prose, triple=True))
        check("the same option in the call IS code",
              "cov_struct=Exchangeable()" in strip_hash_comments(in_code, triple=True))

        print("must_appear")

        def entry(name, expected, files):
            d = Path(tmp) / name
            d.mkdir()
            (d / "expected.json").write_text(json.dumps(expected))
            for fn, body in files.items():
                (d / fn).write_text(body)
            return d

        R_OK = ('# PINNED: ties="efron", because lifelines and coxph must agree on tied event times.\n'
                'fit <- coxph(Surv(time, event) ~ exposed, ties="efron")\n')
        # The SAME file with the option deleted from the statement and left in the comment above
        # it. Byte-identical to R_OK apart from that deletion, which is the point: a substring
        # search over the raw text cannot tell these two apart and calls both a pass.
        R_COMMENT_ONLY = R_OK.replace(', ties="efron")', ')')

        good = entry("good", {"must_appear": {"r": ["coxph(", 'ties="efron"']}},
                     {"r.R": R_OK})
        check("a file containing every declared string passes", lint_must_appear(good) == [],
              lint_must_appear(good))

        comment_only = entry("comment_only", {"must_appear": {"r": ['ties="efron"']}},
                             {"r.R": R_COMMENT_ONLY})
        probs = lint_must_appear(comment_only)
        check("a default named ONLY in a comment is flagged", len(probs) == 1, probs)
        check("and the message says which of the two mistakes it is",
              probs and "only in a comment" in probs[0], probs)

        gone = entry("gone", {"must_appear": {"r": ['ties="efron"']}},
                     {"r.R": "fit <- coxph(Surv(t, e) ~ x)\n"})
        probs = lint_must_appear(gone)
        check("a declared string absent altogether is flagged", len(probs) == 1, probs)
        check("and is reported as absent rather than as a comment",
              probs and "absent" in probs[0], probs)

        # The other direction: a rule that only checks what it was handed asserts its own scope.
        # A file that exists beside the declared one but declares nothing must be flagged.
        undeclared = entry("undeclared", {"must_appear": {"r": ["coxph("]}},
                           {"r.R": "fit <- coxph(Surv(t, e) ~ x)\n",
                            "python.py": "CoxPHFitter().fit(df, 't', 'e')\n"})
        probs = lint_must_appear(undeclared)
        check("a language file that declares NOTHING is flagged", len(probs) == 1, probs)
        check("and the finding names the undeclared language",
              probs and "python" in probs[0], probs)

        ghost = entry("ghost", {"must_appear": {"python": ["CoxPHFitter"]}},
                      {"r.R": "fit <- coxph(Surv(t, e) ~ x)  # no python.py beside it\n"})
        probs = lint_must_appear(ghost)
        check("declaring a language whose file is missing is flagged",
              any("does not exist" in p for p in probs), probs)

        print("public-surface check")
        import public_surface
        # Positive controls: lines exactly as they stood in this repository before the check.
        leaked = [
            "    (n - 1.99) / (n - 1) is exactly 0.99 at n = 100 and the same position at any n (owner decision E6,",
            "when clusters are few — CLAUDE.md warns under 30 and blocks under 15 unless a small-sample approach",
            "This id names the planned page (Ask Mallard docs/specs/code-library-expansion-2026-09-23.md, section 7).",
            "    \"schema prompt warns about: 'silently different defaults are how two correct-looking\",",
            "the one Mallard's evaluation bank calls \"cluster-randomized trial, binary outcome\". Individual-level",
            '"""Controls and fixed-seed calibration for batch F-3.',
            "Claude-Session: https://claude.ai/code/session_01WiV7D9QqVVBRKZ5PUmjPwW",
        ]
        for line in leaked:
            check(f"flags: {line.strip()[:60]}", public_surface.scan_text(line) != [])
        # Negative controls: ordinary lines in this repository that share a word with a pattern.
        clean = [
            '    spec = importlib.util.spec_from_file_location(name, ROOT / "lib" / name / "fixture.py")',
            'print(rbind(sensitivity = sens, specificity = spec))',
            '      "Licensor" shall mean the copyright owner or entity authorized by',
            "    batches = sorted(p for p in here.parent.glob('*_examples.test.py') if p != here)",
            "SEEDS = range(int(os.environ.get(\"MALLARD_SEEDS\", \"50\")))  # 50 in CI; MALLARD_SEEDS=100 reproduces the recorded calibration",
            "# and unweighted tests answer different questions, and a plan that names \"the log-rank test\" means",
        ]
        for line in clean:
            check(f"stays silent: {line.strip()[:60]}", public_surface.scan_text(line) == [])

        print("README list of entries")
        import catalogue
        lib = Path(tmp) / "cat" / "lib"
        for eid, title, eng in [("b-entry", "Beta method", {"r": "executed", "python": "executed"}),
                                ("a-entry", "Alpha method", {"r": "executed", "python": "not-applicable"})]:
            (lib / eid).mkdir(parents=True)
            (lib / eid / "meta.json").write_text(json.dumps({"id": eid, "title": title, "engines": eng}))
        readme = Path(tmp) / "cat" / "README.md"
        readme.write_text(f"intro\n{catalogue.START}\nold\n{catalogue.END}\noutro\n")
        args = ["--readme", str(readme), "--lib", str(lib)]
        check("a stale list fails --check", catalogue.main(args + ["--check"]) == 1)
        catalogue.main(args)
        text = readme.read_text()
        check("after writing, --check passes", catalogue.main(args + ["--check"]) == 0)
        check("text outside the markers is kept", text.startswith("intro\n") and text.endswith("\noutro\n"))
        check("an R-only entry says R only", "| Alpha method | R |" in text, text)
        check("a two-engine entry names both", "| Beta method | R and Python |" in text, text)
        check("the count line counts both kinds",
              "**2 entries.** 1 are executed in both R and Python and carry both claims; 1 are executed in one language" in text, text)
        (lib / "a-entry" / "meta.json").write_text(json.dumps(
            {"id": "a-entry", "title": "Alpha method", "engines": {"r": "executed", "python": "executed"}}))
        check("changing an entry's engines makes the list stale", catalogue.main(args + ["--check"]) == 1)
        # `not-executed` is a state metadata.py accepts (the file exists and CI does not run it). An
        # earlier catalogue rejected it, which would have failed CI on a valid entry.
        (lib / "a-entry" / "meta.json").write_text(json.dumps(
            {"id": "a-entry", "title": "Alpha method", "engines": {"r": "executed", "python": "not-executed"}}))
        catalogue.main(args)
        text = readme.read_text()
        check("a not-executed language is accepted and named as not run",
              "| Alpha method | R (Python present, not run) |" in text, text)
        check("and it does not count as executed", "1 are executed in one language" in text, text)
        (lib / "a-entry" / "meta.json").write_text(json.dumps(
            {"id": "a-entry", "title": "Alpha method", "engines": {"r": "executed", "python": "skipped"}}))
        try:
            catalogue.main(args + ["--check"])
            check("an unknown engine state stops generation", False)
        except SystemExit:
            check("an unknown engine state stops generation", True)

    print()
    if FAILURES:
        print(f"{len(FAILURES)} control(s) failed")
        return 1
    print("all controls passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
