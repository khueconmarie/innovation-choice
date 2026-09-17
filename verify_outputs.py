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
cu=json.loads((P/'checks/characterization_analysis.json').read_text())
assert cu['cases']==48 and cu['max_derivative_error']<1e-8
assert cu['all_signs_checked'] and cu['exact_infinite_geometric_groups']
nu=json.loads((P/'checks/neutrality_analysis.json').read_text())
assert nu['cases']==24 and nu['max_fee_error']<2e-12
assert nu['common_future_invariance'] and nu['min_prefix_investment']>0
assert nu['initial_fee_interiority_checks']==12
assert gu['common_physical_initial_capital_checked']
dy=json.loads((P/'checks/dynamic_analysis.json').read_text())
assert dy['common_physical_initial_capital_checked']
rx=json.loads((P/'checks/resource_example_analysis.json').read_text())
assert rx['cases']==2 and rx['passed'] and rx['precision_decimal_digits']==70
assert rx['max_derivative_error']<1e-20
assert rx['max_double_solver_fee_error']<1e-10
assert rx['same_preferences'] and rx['same_return_paths'] and rx['only_resources_differ']
fg=json.loads((P/'checks/figure_analysis.json').read_text())
assert fg['ordering_checked'] and fg['consumption_rows']==48
for ref in (P/'reference/generated').glob('*.tex'):
    assert ref.read_bytes()==(P/'generated'/ref.name).read_bytes(),ref.name
report=dict(passed=True,csv_rows=csv_rows,npz_arrays=len(a.files),
            max_csv_absolute_difference=max_difference,
            nonlinear_fee_brackets=24,nonlinear_effect_signs=12,
            general_utility_cases=80,characterization_cases=48,neutrality_cases=24,
            compensated_figure_points=48,resource_level_examples=2,initial_fee_interiority_checks=12,
            common_physical_initial_capital_checked=True,
            scope='Output/reference consistency and numerical bound checks; not a new model solution.')
(P/'checks/output_verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
