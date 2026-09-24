"""Controls and fixed-seed calibration for the Bayesian logistic regression entry (R only).

bayesian-logistic-regression. R only and checked against its truth only: Markov chains cannot be
matched draw for draw across languages. The fixture regenerates byte for byte. Oracle: a Laplace
approximation (the posterior mode and the inverse Hessian of the log posterior, computed here with scipy)
must agree with the sampler's posterior mean and SD within Monte Carlo and approximation error. The
sampler's R-hat and effective sample size are checked on the committed fixture. Over the seeds: recovery
of the treatment log odds ratio, and the coverage of the 95% credible interval.

SEEDS: 50 by default, as CI runs it; MALLARD_SEEDS=100 for the recorded calibration. Locally,
MALLARD_R_LIBS may name an R library to put first. Run with --measure to print the numbers.
"""
import csv
import importlib.util
import io
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy import optimize, stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check import parse_harness_block  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
NAME = "bayesian-logistic-regression"
SEEDS = range(int(os.environ.get("MALLARD_SEEDS", "50")))  # 50 in CI; MALLARD_SEEDS=100 reproduces the recorded calibration
PRIOR_SD = np.array([5.0, 2.5, 2.5])


def rank_quantile(misses):
    """The quantile at the same order statistic as the 99th percentile of 100 draws, at any number of draws."""
    n = len(misses)
    return float(np.quantile(misses, (n - 1.99) / (n - 1)))


def load():
    spec = importlib.util.spec_from_file_location(NAME, ROOT / "lib" / NAME / "fixture.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def csv_text(generated):
    out = io.StringIO(newline="")
    writer = csv.DictWriter(out, fieldnames=list(generated[0]))
    writer.writeheader()
    writer.writerows(generated)
    return out.getvalue()


def run_r_many(frames):
    libs = os.environ.get("MALLARD_R_LIBS")
    prelude = f'.libPaths(c("{libs}", .libPaths())); ' if libs else ""
    with tempfile.TemporaryDirectory() as tmp:
        for i, generated in enumerate(frames):
            Path(tmp, str(i)).mkdir()
            Path(tmp, str(i), "fixture.csv").write_text(csv_text(generated))
        script = (ROOT / "lib" / NAME / "r.R").as_posix()
        code = (prelude + f'for (i in 0:{len(frames) - 1}) {{ setwd(file.path("{tmp}", i)); '
                f'out <- capture.output(source("{script}")); cat("@@", i, "\\n"); cat(out, sep = "\\n") }}')
        text = subprocess.check_output(["Rscript", "-e", code], text=True, stderr=subprocess.DEVNULL)
    return [parse_harness_block(chunk) for chunk in text.split("@@")[1:]]


def laplace(generated):
    X = np.array([[1.0, r["treated"], r["age"]] for r in generated])
    y = np.array([r["event"] for r in generated], float)

    def neg(b):
        eta = X @ b
        return -(np.sum(y * eta - np.log1p(np.exp(eta))) - 0.5 * np.sum((b / PRIOR_SD) ** 2))

    mode = optimize.minimize(neg, np.zeros(3), method="BFGS", options={"gtol": 1e-10}).x
    p = 1 / (1 + np.exp(-X @ mode))
    hess = (X * (p * (1 - p))[:, None]).T @ X + np.diag(1 / PRIOR_SD ** 2)
    return mode[1], float(np.sqrt(np.linalg.inv(hess)[1, 1]))


def main():
    m = load()
    committed = list(m.rows())
    assert csv_text(committed).encode() == (ROOT / "lib" / NAME / "fixture.csv").read_bytes(), "fixture.csv does not regenerate"
    results = run_r_many([committed] + [list(m.rows(s)) for s in SEEDS])
    out, runs = results[0], results[1:]
    mode, sd = laplace(committed)
    assert abs(mode - out["treatment_log_or"]) < 0.03, f"posterior mean {out['treatment_log_or']} vs Laplace mode {mode}"
    assert abs(sd / out["posterior_sd"] - 1) < 0.1, f"posterior SD {out['posterior_sd']} vs Laplace {sd}"
    assert out["rhat"] < 1.01 and out["ess"] > 1000, "the sampler has not converged"
    truth = m.truth()["treatment_log_or"]
    coverage = float(np.mean([r["credible_lcl"] <= truth <= r["credible_ucl"] for r in runs]))
    assert coverage >= 0.88, f"credible interval coverage {coverage}"
    assert all(r["rhat"] < 1.02 for r in runs), "a chain failed to mix"
    measured = {"treatment_log_or": rank_quantile([abs(r["treatment_log_or"] - truth) for r in runs]),
                "coverage": coverage, "laplace_mode": float(mode), "laplace_sd": sd}
    if "--measure" in sys.argv:
        print(json.dumps({"seeds": len(SEEDS), "measured": measured}, indent=2))
        return
    exp = json.loads((ROOT / "lib" / NAME / "expected.json").read_text())
    assert abs(exp["truth"]["treatment_log_or"] - truth) < 1e-12
    tol = float(exp["recovery"]["tolerance_estimate"])
    assert measured["treatment_log_or"] < tol, f"recovery miss {measured['treatment_log_or']:.4f} >= {tol}"
    print(json.dumps({"entries": 1, "seeds": len(SEEDS), "result": "pass"}, indent=2))


if __name__ == "__main__":
    main()
