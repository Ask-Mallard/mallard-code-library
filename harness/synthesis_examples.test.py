"""Controls and fixed-seed calibration for the evidence-synthesis entries.

random-effects-meta-analysis, bivariate-diagnostic-meta-analysis (R only), network-meta-analysis. Each
fixture regenerates byte for byte. Oracles: the REML tau^2 maximizes the restricted log-likelihood and
the Hartung-Knapp SE follows its formula; the bivariate model is refitted with metafor::rma.mv (an
independent R implementation of the same REML model); the network estimates are exactly consistent
(C vs B = C vs A - B vs A) and more precise than the direct-only estimate. Over the seeds: the
Hartung-Knapp interval covers the mean effect at least as often as DerSimonian-Laird's z interval; the
network estimator is unbiased. Recovery: each entry on its own seeds (the R-only entry through r.R).

SEEDS: 50 by default, as CI runs it; MALLARD_SEEDS=100 for the recorded calibration. Locally,
MALLARD_R_LIBS may name an R library to put first. Run with --measure to print the numbers.
"""
import contextlib
import csv
import importlib.util
import io
import json
import os
import runpy
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check import parse_harness_block  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
SEEDS = range(int(os.environ.get("MALLARD_SEEDS", "50")))  # 50 in CI; MALLARD_SEEDS=100 reproduces the recorded calibration


def rank_quantile(misses):
    """The quantile at the same order statistic as the 99th percentile of 100 draws, at any number of draws."""
    n = len(misses)
    return float(np.quantile(misses, (n - 1.99) / (n - 1)))


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


def rscript(code):
    libs = os.environ.get("MALLARD_R_LIBS")
    prelude = f'.libPaths(c("{libs}", .libPaths())); ' if libs else ""
    return subprocess.check_output(["Rscript", "-e", prelude + code], text=True, stderr=subprocess.DEVNULL)


def run_r_many(name, frames):
    with tempfile.TemporaryDirectory() as tmp:
        for i, generated in enumerate(frames):
            Path(tmp, str(i)).mkdir()
            Path(tmp, str(i), "fixture.csv").write_text(csv_text(generated))
        script = (ROOT / "lib" / name / "r.R").as_posix()
        text = rscript(f'for (i in 0:{len(frames) - 1}) {{ setwd(file.path("{tmp}", i)); '
                       f'out <- capture.output(source("{script}")); cat("@@", i, "\\n"); cat(out, sep = "\\n") }}')
    return [parse_harness_block(chunk) for chunk in text.split("@@")[1:]]


def recovery(results, truth):
    return {key: rank_quantile([abs(r[key] - value) for r in results]) for key, value in truth.items()}


# ---------------------------------------------------------------------------------------------------

def meta(measured, truths):
    name = "random-effects-meta-analysis"
    regenerates(name)
    m = load(name)
    truths[name] = {"log_or": m.MU}
    out, s = run_python(name)
    f, tau2 = s["neg_restricted_loglik"], out["tau2"]
    assert f(tau2) <= min(f(tau2 * 1.01 + 1e-9), f(max(tau2 * 0.99 - 1e-9, 0))) + 1e-12, "tau^2 is not the REML maximum"
    w, mu = s["pooled"](s["tau2"])
    assert abs(np.sqrt(np.sum(w * (s["y"] - mu) ** 2) / (len(w) - 1) / np.sum(w)) - out["log_or_se"]) < 1e-9
    runs = [run_python(name, list(m.rows(k)))[0] for k in SEEDS]
    z = stats.norm.ppf(0.975)
    hk = np.mean([r["log_or_lcl"] <= m.MU <= r["log_or_ucl"] for r in runs])
    dl = np.mean([abs(r["dl_log_or"] - m.MU) <= z * r["dl_se_z"] for r in runs])
    assert hk >= dl - 0.02 and hk >= 0.88, (hk, dl)
    measured[name] = recovery(runs, truths[name])
    measured[name].update({"coverage_hk": float(hk), "coverage_dl_z": float(dl)})


def dta(measured, truths):
    name = "bivariate-diagnostic-meta-analysis"
    regenerates(name)
    m = load(name)
    truths[name] = m.truth()
    frames = [list(m.rows())] + [list(m.rows(k)) for k in SEEDS]
    results = run_r_many(name, frames)
    committed, runs = results[0], results[1:]
    with tempfile.TemporaryDirectory() as tmp:
        Path(tmp, "fixture.csv").write_text(csv_text(frames[0]))
        text = rscript(f'suppressPackageStartupMessages(library(metafor)); d <- read.csv(file.path("{tmp}", "fixture.csv")); '
                       'long <- rbind(data.frame(study=d$study, o="sens", yi=qlogis(d$TP/(d$TP+d$FN)), vi=1/d$TP+1/d$FN), '
                       'data.frame(study=d$study, o="fpr", yi=qlogis(d$FP/(d$FP+d$TN)), vi=1/d$FP+1/d$TN)); '
                       'f <- rma.mv(yi, vi, mods = ~ o - 1, random = ~ o | study, struct = "UN", data = long, method = "REML", '
                       'control = list(rel.tol = 1e-12)); cat(sprintf("%.12f %.12f", coef(f)[["osens"]], -coef(f)[["ofpr"]]))')
    ls, lp = map(float, text.split())
    assert abs(ls - committed["logit_sensitivity"]) < 1e-5 and abs(lp - committed["logit_specificity"]) < 1e-5, \
        "mada::reitsma and metafor::rma.mv disagree"
    measured[name] = recovery(runs, truths[name])


def nma(measured, truths):
    name = "network-meta-analysis"
    regenerates(name)
    m = load(name)
    truths[name] = m.truth()
    out, s = run_python(name)
    assert abs(out["log_or_C_vs_B"] - (out["log_or_C_vs_A"] - out["log_or_B_vs_A"])) < 1e-9, "estimates must be consistent"
    d = s["d"]
    bc = d[(d.treat1 == "C") & (d.treat2 == "B")]
    w = 1 / bc.seTE.to_numpy() ** 2
    direct_se = float(np.sqrt(1 / w.sum()))
    common_se = float(np.sqrt(np.array([-1.0, 1.0]) @ np.linalg.inv((s["X"] * s["w"][:, None]).T @ s["X"]) @ np.array([-1.0, 1.0])))
    assert common_se < direct_se, "indirect evidence should add precision"
    runs = [run_python(name, list(m.rows(k)))[0] for k in SEEDS]
    for key, value in truths[name].items():
        assert abs(np.mean([r[key] for r in runs]) - value) < 0.05, f"{key} biased"
    measured[name] = recovery(runs, truths[name])
    measured[name].update({"direct_se_C_vs_B": direct_se, "network_common_se_C_vs_B": common_se})


def check_tolerances(measured, truths):
    for name, truth in truths.items():
        exp = json.loads((ROOT / "lib" / name / "expected.json").read_text())
        for key, value in truth.items():
            assert abs(exp["truth"][key] - value) < 1e-10, f"{name} {key}: expected.json {exp['truth'][key]} != {value}"
            tol = float(exp["recovery"].get("tolerance_by_key", {}).get(key, exp["recovery"]["tolerance_estimate"]))
            assert measured[name][key] < tol, f"{name} {key}: recovery miss {measured[name][key]:.4f} >= {tol}"


if __name__ == "__main__":
    measured, truths = {}, {}
    meta(measured, truths)
    dta(measured, truths)
    nma(measured, truths)
    if "--measure" in sys.argv:
        print(json.dumps({"seeds": len(SEEDS), "measured": measured, "truths": truths}, indent=2))
        sys.exit(0)
    check_tolerances(measured, truths)
    print(json.dumps({"entries": len(measured), "seeds": len(SEEDS), "result": "pass"}, indent=2))
