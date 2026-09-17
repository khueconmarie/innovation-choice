"""Regenerate all numerical results, then compare them with release references."""
from pathlib import Path
import json, os, platform, subprocess, sys, time
from importlib.metadata import version

P = Path(__file__).resolve().parent
for d in ['checks','generated','figures','logs']:
    (P/d).mkdir(exist_ok=True)
scripts = ['general_utility_analysis.py','characterization_analysis.py','neutrality_analysis.py','resource_example.py',
           'dynamic_analysis.py','scale_analysis.py',
           'joint_analysis.py','nonlinear_analysis.py','render_figure.py','verify_outputs.py']
record = dict(python=sys.version, platform=platform.platform(),
              packages={p:version(p) for p in ['numpy','scipy','sympy','matplotlib','mpmath']},
              scripts=[])
for script in scripts:
    start = time.monotonic()
    print(f'Running {script}', flush=True)
    with (P/'logs'/f'{script}.log').open('w') as out:
        result = subprocess.run([sys.executable, str(P/script)], cwd=P,
                                stdout=out, stderr=subprocess.STDOUT,
                                env=dict(os.environ, MPLBACKEND='Agg'))
    record['scripts'].append(dict(script=script, exit_code=result.returncode,
                                 seconds=time.monotonic()-start))
    if result.returncode:
        print((P/'logs'/f'{script}.log').read_text())
        raise SystemExit(result.returncode)
record['total_seconds'] = sum(s['seconds'] for s in record['scripts'])
record['passed'] = True
(P/'checks/reproduction_run.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
