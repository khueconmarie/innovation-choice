"""Exact infinite-stream checks of the universal-curvature necessity witness.

Only three consumption groups are needed after a single common-return shock.
No terminal approximation is used. Numerical checks supplement the proof.
"""
from pathlib import Path
import csv, json
import numpy as np
from scipy.optimize import brentq
from general_utility_analysis import Utility
P=Path(__file__).resolve().parent


def allocation(u, beta, gain, shock, *, resources=None, target=None):
    # Relative to R_t=1/beta, group factors are 1, gain, gain*exp(shock).
    factors=np.array([1.,gain,gain*np.exp(shock)])
    utility_weights=np.array([1.,beta,beta**2/(1-beta)])
    prices=utility_weights/factors
    def values(ll):
        c=u.demand(ll-np.log(factors))
        return c, float(prices@c), float(utility_weights@u.value(c))
    index,goal=(1,resources) if resources is not None else (2,target)
    ll=brentq(lambda x: values(x)[index]-goal,-40.,40.,xtol=1e-14)
    c,pv,v=values(ll)
    return dict(c=c,pv=pv,value=v,lam=np.exp(ll))


def pair(u,beta,c,eps,shock=0.):
    X=c/(1-beta)
    outside=allocation(u,beta,1.,shock,resources=X)
    access=allocation(u,beta,1+eps,shock,target=outside['value'])
    return X-access['pv'],outside,access


def main():
    rows=[]
    for name,gs,ws in [('low',[.3,.7],[.4,.6]),('high',[1.3,1.7],[.6,.4]),
                       ('log',[1.],[1.]),('crossing',[.5,1.5],[.5,.5])]:
        u=Utility(gs,ws)
        for beta in [.9,.96]:
            for c in [.4,1.4,3.]:
                curvature=sum(g*w*c**(-g) for g,w in zip(gs,ws))/u.marginal(c)
                for eps in [.001,.01]:
                    fee,b,a=pair(u,beta,c,eps)
                    assert 0<fee<c/(1-beta)
                    assert a['c'][0]<c and np.all(a['c'][1:]>c)
                    derivative=beta**2/(1-beta)/a['lam']*(a['c'][2]*u.marginal(a['c'][2])-c*u.marginal(c))
                    h=2e-4
                    fd=(pair(u,beta,c,eps,h)[0]-pair(u,beta,c,eps,-h)[0])/(2*h)
                    if name=='log': assert abs(derivative)<1e-12
                    else: assert derivative*(1-curvature)>0
                    rows.append(dict(utility=name,beta=beta,c=c,epsilon=eps,curvature=float(curvature),
                        fee=fee,c_initial=float(a['c'][0]),c_later=float(a['c'][2]),
                        derivative=float(derivative),finite_difference=float(fd),
                        derivative_error=float(abs(derivative-fd)),
                        compensation_error=abs(a['value']-b['value'])))
    err=max(x['derivative_error'] for x in rows)
    assert err<1e-8,err
    (P/'checks').mkdir(exist_ok=True);(P/'generated').mkdir(exist_ok=True)
    with (P/'generated/characterization_witnesses.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    report=dict(cases=len(rows),max_derivative_error=err,
        max_compensation_error=max(x['compensation_error'] for x in rows),
        all_signs_checked=True,exact_infinite_geometric_groups=True,
        scope='Numerical checks of constructive necessity; the proof establishes universal quantifiers.')
    (P/'checks/characterization_analysis.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__': main()
