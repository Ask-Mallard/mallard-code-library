"""Controls and fixed-seed calibration for batch F-2.

difference-in-differences-twfe, interrupted-time-series-segmented, controlled-interrupted-time-series,
regression-discontinuity-local-linear. Each fixture regenerates byte for byte. Exact oracles: the 2x2 DiD
of group means equals the two-way fixed-effects coefficient; the Newey-West covariance is rebuilt by
hand; the difference-series coefficients equal the stacked interaction model's; the conventional RD
estimate equals a hand-written triangular-kernel local linear fit at rdrobust's bandwidth. Each entry's
naive comparator is shown to miss. Recovery is calibrated by running each entry's own python.py on 100
generated fixtures. Run with --measure to print the numbers.
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
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check import parse_harness_block  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
SEEDS = range(100)


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "lib" / name / "fixture.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def csv_text(generated):
    out = io.StringIO(newline="")
    writer = csv.DictWriter(out, fieldnames=list(generated[0]))
    writer.writeheader()
    writer.writerows(generated)
    return out.getvalue()


def regenerates(name):
    assert csv_text(list(load(name).rows())).encode() == (ROOT / "lib" / name / "fixture.csv").read_bytes(), \
        f"{name}: fixture.csv does not regenerate byte for byte"


def run_python(name, generated=None):
    script = ROOT / "lib" / name / "python.py"
    with tempfile.TemporaryDirectory() as tmp:
        src = (ROOT / "lib" / name / "fixture.csv").read_text() if generated is None else csv_text(generated)
        Path(tmp, "fixture.csv").write_text(src)
        here = os.getcwd()
        buf = io.StringIO()
        try:
            os.chdir(tmp)
            with contextlib.redirect_stdout(buf):
                scope = runpy.run_path(str(script))
        finally:
            os.chdir(here)
    return parse_harness_block(buf.getvalue()), scope


def calibrate(name, truth):
    module = load(name)
    runs = [run_python(name, list(module.rows(s)))[0] for s in SEEDS]
    return {key: float(np.quantile([abs(r[key] - value) for r in runs], 0.99)) for key, value in truth.items()}


def segmented(d, ycol):
    X = np.column_stack([np.ones(len(d)), d.month, d.post, d.months_since]).astype(float)
    return X, d[ycol].to_numpy(float)


def newey_west(X, resid, lag):
    S = (X * resid[:, None]).T @ (X * resid[:, None])
    for j in range(1, lag + 1):
        G = (X[j:] * resid[j:, None]).T @ (X[:-j] * resid[:-j, None])
        S += (1 - j / (lag + 1)) * (G + G.T)
    bread = np.linalg.inv(X.T @ X)
    return bread @ S @ bread


# ---------------------------------------------------------------------------------------------------

def did(measured, truths):
    name = "difference-in-differences-twfe"
    regenerates(name)
    m = load(name)
    truths[name] = {"did": m.EFFECT}
    out, scope = run_python(name)
    d = scope["d"]
    g = lambda a, p: d.loc[(d.adopter == a) & (d.post == p), "los"].mean()
    two_by_two = (g(1, 1) - g(1, 0)) - (g(0, 1) - g(0, 0))
    assert abs(two_by_two - out["did"]) < 1e-10, "TWFE is not the 2x2 DiD on a balanced single-date panel"
    assert abs(out["pre_post_adopters"] - m.EFFECT) > 1.0 and abs(out["post_between_groups"] - m.EFFECT) > 1.0
    measured[name] = calibrate(name, truths[name])


def its(measured, truths):
    name = "interrupted-time-series-segmented"
    regenerates(name)
    m = load(name)
    truths[name] = {"level_change": m.LEVEL_CHANGE, "slope_change": m.SLOPE_CHANGE}
    out, scope = run_python(name)
    X, y = segmented(scope["d"], "rate")
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    V = newey_west(X, y - X @ beta, 3)
    assert abs(np.sqrt(V[2, 2]) - out["level_change_se"]) < 1e-10
    assert abs(np.sqrt(V[3, 3]) - out["slope_change_se"]) < 1e-10
    # Positive autocorrelation: the OLS standard errors are the smaller ones.
    assert out["residual_lag1_autocorrelation"] > 0.2
    assert out["level_change_se_ols"] < 0.8 * out["level_change_se"]
    assert out["slope_change_se_ols"] < 0.8 * out["slope_change_se"]
    measured[name] = calibrate(name, truths[name])


def cits(measured, truths):
    name = "controlled-interrupted-time-series"
    regenerates(name)
    m = load(name)
    truths[name] = {"level_change": m.LEVEL_CHANGE, "slope_change": m.SLOPE_CHANGE}
    out, scope = run_python(name)
    d = scope["d"]
    # The stacked interaction model: its group-by-period coefficients equal the difference series'.
    stacked = pd.concat([
        d.assign(y=d.intervention_rate, g=1.0), d.assign(y=d.control_rate, g=0.0)], ignore_index=True)
    X = np.column_stack([np.ones(len(stacked)), stacked.month, stacked.post, stacked.months_since,
                         stacked.g, stacked.g * stacked.month, stacked.g * stacked.post,
                         stacked.g * stacked.months_since]).astype(float)
    beta = np.linalg.lstsq(X, stacked.y.to_numpy(float), rcond=None)[0]
    assert abs(beta[6] - out["level_change"]) < 1e-9 and abs(beta[7] - out["slope_change"]) < 1e-9
    Xd, yd = segmented(d, "difference")
    bd = np.linalg.lstsq(Xd, yd, rcond=None)[0]
    assert abs(np.sqrt(newey_west(Xd, yd - Xd @ bd, 3)[2, 2]) - out["level_change_se"]) < 1e-10
    # The single-series analysis credits the shared drop to the intervention.
    assert out["single_series_level_change"] - m.LEVEL_CHANGE < -1.0
    measured[name] = calibrate(name, truths[name])


def rd(measured, truths):
    name = "regression-discontinuity-local-linear"
    regenerates(name)
    m = load(name)
    truths[name] = {"jump": m.JUMP}
    out, scope = run_python(name)
    d = scope["d"]
    h = out["bandwidth_h"]
    x, y = d.score.to_numpy(float), d.y.to_numpy(float)

    def side(mask):
        w = np.clip(1 - np.abs(x[mask]) / h, 0, None)
        keep = w > 0
        X = np.column_stack([np.ones(keep.sum()), x[mask][keep]])
        W = w[keep]
        return np.linalg.solve(X.T @ (X * W[:, None]), X.T @ (W * y[mask][keep]))[0], keep.sum()

    right, n_right = side(x >= 0)
    left, n_left = side(x < 0)
    assert abs((right - left) - out["jump"]) < 1e-8, "conventional estimate is not the local linear fit"
    assert (n_left, n_right) == (out["n_left_in_h"], out["n_right_in_h"])
    assert abs(out["global_linear_jump"] - m.JUMP) > 0.5
    assert out["robust_ucl"] - out["robust_lcl"] > out["conventional_ucl"] - out["conventional_lcl"]
    measured[name] = calibrate(name, truths[name])


def check_tolerances(measured, truths):
    for name, truth in truths.items():
        exp = json.loads((ROOT / "lib" / name / "expected.json").read_text())
        for key, value in truth.items():
            assert abs(exp["truth"][key] - value) < 1e-12, f"{name} {key}: expected.json {exp['truth'][key]} != {value}"
            tol = float(exp["recovery"].get("tolerance_by_key", {}).get(key, exp["recovery"]["tolerance_estimate"]))
            assert measured[name][key] < tol, f"{name} {key}: 99th percentile miss {measured[name][key]:.4f} >= {tol}"


if __name__ == "__main__":
    measured, truths = {}, {}
    did(measured, truths)
    its(measured, truths)
    cits(measured, truths)
    rd(measured, truths)
    if "--measure" in sys.argv:
        print(json.dumps({"measured": measured, "truths": truths}, indent=2))
        sys.exit(0)
    check_tolerances(measured, truths)
    print(json.dumps({"entries": len(measured), "seeds": len(SEEDS), "result": "pass"}, indent=2))
