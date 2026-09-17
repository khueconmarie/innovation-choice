"""Exact infinite-stream resource comparison with independent high precision.

Same u(c)=sqrt(c)-1/sqrt(c), beta, and return paths in both rows. Only X
changes. Arbitrary-precision compensation equations and numerical derivatives
are checked against the separately implemented three-group allocation solver.
"""
from pathlib import Path
import csv,json
import mpmath as mp
from characterization_analysis import pair
from general_utility_analysis import Utility

P=Path(__file__).resolve().parent
mp.mp.dps=70
b=mp.mpf('.96');gain=mp.mpf('1.01')
u=lambda c:mp.sqrt(c)-1/mp.sqrt(c)
du=lambda c:(c+1)/(2*c**mp.mpf('1.5'))
weights=[mp.mpf(1),b,b*b/(1-b)]

def solve(c,shock=mp.mpf(0)):
    X=c/(1-b)
    factors=[mp.mpf(1),mp.mpf(1),mp.exp(shock)]
    def outside_eq(x0,x1,x2):
        xs=[x0,x1,x2]
        return (du(x0)-factors[1]*du(x1),du(x0)-factors[2]*du(x2),
                sum(w*x/f for w,x,f in zip(weights,xs,factors))-X)
    cb=mp.findroot(outside_eq,(c,c,c),tol=mp.mpf('1e-60'))
    target=sum(w*u(x) for w,x in zip(weights,cb))
    fi=[1,gain,gain*mp.exp(shock)]
    def access_eq(x0,x1,x2):
        return (du(x0)-fi[1]*du(x1),du(x0)-fi[2]*du(x2),
                sum(w*u(x) for w,x in zip(weights,[x0,x1,x2]))-target)
    ca=mp.findroot(access_eq,(c*.99,c*1.001,c*1.001),tol=mp.mpf('1e-60'))
    cost=sum(w*x/f for w,x,f in zip(weights,ca,fi))
    return X-cost,cb,ca

rows=[]
for cs in ['.4','3']:
    c=mp.mpf(cs);X=c/(1-b);fee,cb,ca=solve(c)
    derivative=b*b/((1-b)*du(ca[0]))*(ca[2]*du(ca[2])-c*du(c))/X
    h=mp.mpf('1e-12')
    numerical=(solve(c,h)[0]-solve(c,-h)[0])/(2*h*X)
    assert abs(derivative-numerical)<mp.mpf('1e-25')
    f2,_,a2=pair(Utility([.5,1.5],[.5,.5]),float(b),float(c),.01)
    assert abs(float(fee)-f2)<2e-12
    assert ca[0]<c<ca[2]
    assert (ca[2]<1) if c<1 else (c>1)
    rows.append(dict(c=float(c),X=float(X),beta=float(b),epsilon=.01,
                     curvature=float((c+3)/(2*(c+1))),fee_fraction=float(fee/X),
                     c_initial=float(ca[0]),c_later=float(ca[2]),
                     derivative_fee_fraction=float(derivative),
                     finite_difference=float(numerical),
                     high_precision_derivative=mp.nstr(derivative,45),
                     finite_difference_error=float(abs(derivative-numerical)),
                     double_solver_fee_error=abs(float(fee)-f2)))
assert rows[0]['derivative_fee_fraction']<0<rows[1]['derivative_fee_fraction']
for name in ['generated','checks']:(P/name).mkdir(exist_ok=True)
with (P/'generated/resource_level_example.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
report=dict(passed=True,cases=2,precision_decimal_digits=70,same_preferences=True,
            same_return_paths=True,only_resources_differ=True,exact_infinite_groups=True,
            max_derivative_error=max(r['finite_difference_error'] for r in rows),
            max_double_solver_fee_error=max(r['double_solver_fee_error'] for r in rows),
            scope='Two local theoretical examples; no empirical wealth threshold or global monotone sign-switch claim.')
(P/'checks/resource_example_analysis.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'checks':report,'rows':rows},indent=2))
