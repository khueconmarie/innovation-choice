"""Compare generated output with reference data; this does not solve a model."""
from pathlib import Path
import csv, json
import numpy as np

P = Path(__file__).resolve().parent
errors = []
max_difference = 0.
csv_rows = 0
for ref in (P/'reference/generated').glob('*.csv'):
    actual = P/'generated'/ref.name
    with ref.open() as f: a = list(csv.DictReader(f))
    with actual.open() as f: b = list(csv.DictReader(f))
    assert len(a) == len(b),ref.name
    csv_rows += len(a)
    for x,y in zip(a,b):
        assert set(x)==set(y)
        for key in x:
            try:
                n,m=float(x[key]),float(y[key])
            except ValueError:
                assert x[key]==y[key],(ref.name,key)
            else:
                assert np.isclose(n,m,rtol=2e-8,atol=2e-10),(ref.name,key,n,m)
                max_difference=max(max_difference,abs(n-m))
a=np.load(P/'reference/generated/nonlinear_paths.npz')
b=np.load(P/'generated/nonlinear_paths.npz')
assert set(a.files)==set(b.files)
for key in a.files:
    assert np.allclose(a[key],b[key],rtol=2e-8,atol=2e-10),key
nl=json.loads((P/'checks/nonlinear_analysis.json').read_text())
assert len(nl['cases'])==24 and len(nl['extended_horizon'])==24
assert all(x['change_interval'][0]*x['change_interval'][1]>0 for x in nl['sign_intervals'])
assert nl['max_fee_interval_width']<3e-8
assert nl['max_euler_log_residual']<5e-11
assert nl['max_horizon_fee_change']<3e-10
assert nl['max_envelope_error']<5e-8
gu=json.loads((P/'checks/general_utility_analysis.json').read_text())
assert gu['cases']==80 and gu['max_derivative_error']<3e-8
assert gu['min_resource_cap_slack']>0
report=dict(passed=True,csv_rows=csv_rows,npz_arrays=len(a.files),
            max_csv_absolute_difference=max_difference,
            nonlinear_fee_brackets=24,nonlinear_effect_signs=12,
            general_utility_cases=80,
            scope='Output/reference consistency and numerical bound checks; not a new model solution.')
(P/'checks/output_verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
