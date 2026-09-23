"""Numerical positive/negative controls and fixed-seed calibration for the new examples.

The small survey oracle uses exact rational arithmetic, separate from either engine.
AUC's oracle counts every case-control pair, separate from the rank implementations.
"""
import importlib.util
import json
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.special import logit
import statsmodels.api as sm

ROOT = Path(__file__).resolve().parent.parent


def fixture(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'lib' / name / 'fixture.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.rows


def survey(rows, domain=False, broken=False):
    if broken:
        rows = [r for r in rows if r['domain']]
    total = sum(F(r['weight']) for r in rows if not domain or r['domain'])
    mean = sum(F(r['weight'])*r['outcome'] for r in rows if not domain or r['domain'])/total
    psus = {}
    for r in rows:
        key = (r['stratum'], r['psu'])
        psus[key] = psus.get(key, F(0)) + F(r['weight'])*(r['domain'] if domain else 1)*(r['outcome']-mean)/total
    variance = F(0)
    for h in {h for h, _ in psus}:
        values = [v for (hh, _), v in psus.items() if hh == h]
        m = len(values)
        assert m > 1
        center = sum(values)/m
        variance += F(m, m-1)*sum((v-center)**2 for v in values)
    return float(mean), float(variance)**.5


if __name__ == '__main__':
    import csv
    import io
    from check import parse_harness_block
    import subprocess
    # Regeneration must reproduce each committed fixture byte for byte.
    for name in ('complex-survey-domain-prevalence', 'prediction-binary-validation'):
        generated = list(fixture(name)())
        output = io.StringIO(newline='')
        writer = csv.DictWriter(output, fieldnames=list(generated[0]))
        writer.writeheader()
        writer.writerows(generated)
        assert output.getvalue().encode() == (ROOT/'lib'/name/'fixture.csv').read_bytes(), name
    generate = fixture('complex-survey-domain-prevalence')
    rows = list(generate())
    full, domain = survey(rows), survey(rows, True)
    wrong = survey(rows, True, True)
    assert abs(domain[1]-wrong[1]) > 1e-6, 'dropping empty-domain PSUs must be detectable'
    # Exact byte-faithful committed CSV rather than an approximation of the negative fixture.
    path = ROOT / 'lib/complex-survey-domain-prevalence'
    actual = parse_harness_block(subprocess.check_output([__import__('sys').executable, 'python.py'], cwd=path, text=True))
    assert abs(actual['prevalence']-full[0]) < 1e-9
    assert abs(actual['domain_prevalence_se']-domain[1]) < 1e-9
    truth = (.15+2*.35+4*.60)/7
    misses = [max(abs(survey(list(generate(seed)))[0]-truth), abs(survey(list(generate(seed)), True)[0]-truth)) for seed in range(100)]
    assert np.quantile(misses,.99) < .08
    # Raw validation observations for an independent pairwise AUC oracle.
    path = ROOT / 'lib/prediction-binary-validation'
    records = list(csv.DictReader((path/'fixture.csv').open()))
    p=np.array([float(r['predicted']) for r in records]); y=np.array([int(r['outcome']) for r in records])
    pairs = p[y==1,None] - p[y==0][None,:]
    auc = float(np.mean((pairs>0)+.5*(pairs==0)))
    actual = parse_harness_block(subprocess.check_output([__import__('sys').executable,'python.py'],cwd=path,text=True))
    assert abs(actual['auc']-auc) < 1e-9
    assert abs(actual['auc']-(1-auc)) > .1, 'reversing event direction must be detectable'
    assert abs(actual['brier']-np.mean((y-p)**2)) < 1e-9
    generate = fixture('prediction-binary-validation')
    calibration = []
    for seed in range(100):
        rows=list(generate(seed));p=np.array([r['predicted'] for r in rows]);y=np.array([r['outcome'] for r in rows]);lp=logit(p)
        fit=sm.GLM(y,sm.add_constant(lp),family=sm.families.Binomial()).fit()
        assert fit.converged
        calibration.append(max(abs(fit.params[0]),abs(fit.params[1]-1)))
    assert np.quantile(calibration,.99) < .25
    print(json.dumps({'survey':{'se_correct':domain[1],'se_deleted_psus':wrong[1],'recovery_miss_99th':float(np.quantile(misses,.99))},'prediction':{'pairwise_auc':auc,'calibration_miss_99th':float(np.quantile(calibration,.99))},'seeds':100,'result':'pass'},indent=2))
    # The descriptive entries' controls run from here because CI's controls step already calls this
    # file; listing descriptive_examples.test.py in verify.yml directly needs a workflow-scoped
    # change, and a control CI never runs is decoration. check_call fails this file if they fail.
    subprocess.check_call([__import__('sys').executable, str(Path(__file__).with_name('descriptive_examples.test.py'))])
