"""Controls and fixed-seed calibration for batch H-3 (penalized development with optimism correction).

penalized-prediction-optimism. The fixture (fixture.csv and bootstrap.csv) regenerates byte for byte.
Oracles: the ridge fit satisfies its penalized score equations; with a vanishing penalty it equals the
unpenalized logistic MLE (statsmodels); the selected penalty is the argmin of the cross-validated
deviance. Over the seeds: the apparent AUC exceeds the external AUC on average (optimism), and the
bootstrap-corrected AUC is closer to it. Recovery: the corrected-minus-external AUC, whose truth is 0.

SEEDS: 50 by default, as CI runs it; MALLARD_SEEDS=100 for the recorded calibration (owner decisions E5,
E6). Run with --measure to print the numbers.
"""
import contextlib
import csv
import importlib.util
import io
import json
import os
import runpy
import sys
import tempfile
from pathlib import Path

import numpy as np
import statsmodels.api as sm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check import parse_harness_block  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
NAME = "penalized-prediction-optimism"
SEEDS = range(int(os.environ.get("MALLARD_SEEDS", "50")))  # 50 in CI; MALLARD_SEEDS=100 reproduces the recorded calibration (owner decisions E5, E6)


def rank_quantile(misses):
    """The quantile at the same order statistic as the 99th percentile of 100 draws (owner decision E6)."""
    n = len(misses)
    return float(np.quantile(misses, (n - 1.99) / (n - 1)))


def load():
    spec = importlib.util.spec_from_file_location(NAME, ROOT / "lib" / NAME / "fixture.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def texts(m, seed):
    generated = list(m.rows(seed))
    out = io.StringIO(newline="")
    writer = csv.DictWriter(out, fieldnames=list(generated[0]))
    writer.writeheader()
    writer.writerows(generated)
    boot = io.StringIO(newline="")
    w = csv.writer(boot)
    w.writerow(["resample", "rows"])
    for b, idx in enumerate(m.bootstrap(seed), start=1):
        w.writerow([b, " ".join(map(str, idx))])
    return out.getvalue(), boot.getvalue()


def run_python(fixture_text, boot_text):
    with tempfile.TemporaryDirectory() as tmp:
        Path(tmp, "fixture.csv").write_text(fixture_text)
        Path(tmp, "bootstrap.csv").write_text(boot_text)
        here = os.getcwd()
        buf = io.StringIO()
        try:
            os.chdir(tmp)
            with contextlib.redirect_stdout(buf):
                scope = runpy.run_path(str(ROOT / "lib" / NAME / "python.py"))
        finally:
            os.chdir(here)
    return parse_harness_block(buf.getvalue()), scope


def main():
    m = load()
    fixture_text, boot_text = texts(m, m.SEED)
    assert fixture_text.encode() == (ROOT / "lib" / NAME / "fixture.csv").read_bytes(), "fixture.csv does not regenerate"
    assert boot_text.encode() == (ROOT / "lib" / NAME / "bootstrap.csv").read_bytes(), "bootstrap.csv does not regenerate"
    out, s = run_python(fixture_text, boot_text)

    # The ridge fit solves its penalized score equations on the standardized scale.
    X, y, lam, fit = s["X"], s["y"], s["lam"], s["fit"]
    mu, sd = X.mean(axis=0), X.std(axis=0, ddof=1)
    Z = np.column_stack([np.ones(len(X)), (X - mu) / sd])
    b_std = np.r_[fit[0] + np.sum(fit[1:] * mu), fit[1:] * sd]
    p = 1 / (1 + np.exp(-Z @ b_std))
    grad = Z.T @ (y - p) - np.r_[0.0, lam * b_std[1:]]
    assert np.max(np.abs(grad)) < 1e-8, "the ridge fit does not solve its score equations"
    # A vanishing penalty gives the unpenalized maximum likelihood fit.
    mle = sm.Logit(y, sm.add_constant(X)).fit(disp=0, tol=1e-12).params
    assert np.max(np.abs(s["ridge"](X, y, 1e-10) - mle)) < 1e-6
    # The selected penalty minimizes the cross-validated deviance over the grid.
    fold = s["dev"]["fold"].to_numpy()
    cv = [sum(s["deviance"](s["lp"](s["ridge"](X[fold != f], y[fold != f], l), X[fold == f]), y[fold == f])
              for f in np.unique(fold)) for l in s["lambdas"]]
    assert s["lambdas"][int(np.argmin(cv))] == lam

    runs = [run_python(*texts(m, seed))[0] for seed in SEEDS]
    app = np.mean([r["apparent_auc"] - r["external_auc"] for r in runs])
    cor = np.mean([r["corrected_auc"] - r["external_auc"] for r in runs])
    assert app > 0.02, "apparent AUC should be optimistic"
    assert abs(cor) < abs(app) / 2, "the corrected AUC should be much closer to the external one"
    app_s = np.mean([r["apparent_slope"] - r["external_slope"] for r in runs])
    cor_s = np.mean([r["corrected_minus_external_slope"] for r in runs])
    assert app_s > 0.1 and abs(cor_s) < abs(app_s) / 2, "the corrected slope should remove the optimism"
    measured = {"corrected_minus_external_auc": rank_quantile([abs(r["corrected_minus_external_auc"]) for r in runs]),
                "mean_apparent_minus_external_auc": float(app), "mean_corrected_minus_external_auc": float(cor),
                "mean_apparent_minus_external_slope": float(np.mean([r["apparent_slope"] - r["external_slope"] for r in runs])),
                "mean_corrected_minus_external_slope": float(np.mean([r["corrected_minus_external_slope"] for r in runs])),
                "corrected_minus_external_slope": rank_quantile([abs(r["corrected_minus_external_slope"]) for r in runs])}
    if "--measure" in sys.argv:
        print(json.dumps({"seeds": len(SEEDS), "measured": measured}, indent=2))
        return
    exp = json.loads((ROOT / "lib" / NAME / "expected.json").read_text())
    for key in ("corrected_minus_external_auc", "corrected_minus_external_slope"):
        tol = float(exp["recovery"]["tolerance_by_key"][key])
        assert exp["truth"][key] == 0.0
        assert measured[key] < tol, f"{key}: recovery miss {measured[key]:.4f} >= {tol}"
    print(json.dumps({"entries": 1, "seeds": len(SEEDS), "result": "pass"}, indent=2))


if __name__ == "__main__":
    main()
