"""Independent policy/value checks of two logarithmic-neutrality formulas.

Simulates a finite nonstationary prefix and appends its exact stationary
infinite continuation; the fee is solved from the resulting welfare equality.
"""
from pathlib import Path
import csv, json
import numpy as np
from scipy.optimize import brentq
P=Path(__file__).resolve().parent


def value_cobb(beta,eta,A,K,qs,tail_q,fee=0.):
    b=eta/(1-beta*eta)
    def constant(q):
        return (np.log(1-beta*eta)+(1+beta*b)*np.log(A)+beta*b*np.log(beta*eta*q))/(1-beta)
    v=0.; minj=float('inf')
    for t,q in enumerate(qs):
        resources=A*K**eta-(fee if t==0 else 0.)
        if resources<=0: return -np.inf,0.
        J=beta*eta*resources; C=resources-J
        v+=beta**t*np.log(C); K=q*J; minj=min(minj,J)
    v+=beta**len(qs)*(b*np.log(K)+constant(tail_q))
    return float(v),minj


def value_linear_survival(beta,sigma,K,qs,tail_q,fee=0.):
    a=1/(1-beta)
    d=lambda q: np.log(1-beta)+beta*a*np.log(beta)+a*np.log(sigma+q)-np.log(q)
    v=0.; minj=float('inf')
    for t,q in enumerate(qs):
        f=fee if t==0 else 0.
        W=(1+sigma/q)*K-f
        if W<=0: return -np.inf,0.
        C=(1-beta)*W;J=K-f-C
        assert J>0,('Interior domain failed',beta,sigma,t,fee,J)
        v+=beta**t*np.log(C);K=sigma*K+q*J;minj=min(minj,J)
    assert beta*tail_q>(1-beta)*sigma
    v+=beta**len(qs)*(a*np.log(K)+d(tail_q)/(1-beta))
    return float(v),minj


def main():
    rng=np.random.default_rng(420917)
    rows=[];T=5;N=12;K0=1.3
    for beta in [.90,.96]:
        base=rng.uniform(1.005,1.045,T)
        access=base*np.exp(rng.uniform(.01,.03,T))
        for eta in [.36,.9,1.]:
            A=.8;Y0=A*K0**eta
            formula=Y0*(-np.expm1(-eta*np.sum(beta**(np.arange(T)+1)*np.log(access/base))))
            for common in [1.05,1.12]:
                future=np.linspace(common*.98,common,N-T)
                qb=np.r_[base,future];qi=np.r_[access,future]
                vb,_=value_cobb(beta,eta,A,K0,qb,common)
                numerical=brentq(lambda f:value_cobb(beta,eta,A,K0,qi,common,f)[0]-vb,0.,Y0*.4,xtol=1e-14)
                vi,minj=value_cobb(beta,eta,A,K0,qi,common,numerical)
                rows.append(dict(case='cobb_full_depreciation',beta=beta,eta=eta,sigma=0.,common_q=common,
                    formula_fee=formula,numerical_fee=numerical,fee_error=abs(formula-numerical),
                    value_error=abs(vi-vb),min_prefix_investment=minj))
        for sigma in [.0,.2,.9]:
            a=1/(1-beta)
            deltaB=np.sum(beta**np.arange(T)*(a*np.log((sigma+access)/(sigma+base))-np.log(access/base)))
            formula=K0*(1+sigma/access[0])*(-np.expm1(-(1-beta)*deltaB))
            for common in [1.05,1.12]:
                future=np.linspace(common*.98,common,N-T)
                qb=np.r_[base,future];qi=np.r_[access,future]
                vb,_=value_linear_survival(beta,sigma,K0,qb,common)
                numerical=brentq(lambda f:value_linear_survival(beta,sigma,K0,qi,common,f)[0]-vb,0.,K0*.4,xtol=1e-14)
                vi,minj=value_linear_survival(beta,sigma,K0,qi,common,numerical)
                rows.append(dict(case='linear_with_survival',beta=beta,eta=1.,sigma=sigma,common_q=common,
                    formula_fee=formula,numerical_fee=numerical,fee_error=abs(formula-numerical),
                    value_error=abs(vi-vb),min_prefix_investment=minj))
    for i in range(0,len(rows),2):
        assert abs(rows[i]['numerical_fee']-rows[i+1]['numerical_fee'])<2e-12
    error=max(r['fee_error'] for r in rows);assert error<2e-12,error
    (P/'generated').mkdir(exist_ok=True);(P/'checks').mkdir(exist_ok=True)
    with (P/'generated/neutrality_boundaries.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    report=dict(cases=len(rows),max_fee_error=error,
        max_value_error=max(r['value_error'] for r in rows),
        min_prefix_investment=min(r['min_prefix_investment'] for r in rows),
        common_future_invariance=True,exact_tail=True,
        scope='Independent simulation and exact infinite continuation under the maintained interior domain.')
    # Directly verify the initial-fee condition separately from no-fee dates.
    margins=[]
    for beta in [.90,.96]:
        for sigma in [.0,.2,.9]:
            q=1.03;K=1.3
            limit=1-(1-beta)*sigma/(beta*q)
            for fraction in [.25,.75]:
                fee=K*limit*fraction
                W=(1+sigma/q)*K-fee
                direct=K-fee-(1-beta)*W
                formula=(beta-(1-beta)*sigma/q)*K-beta*fee
                assert abs(direct-formula)<1e-14 and direct>0
                assert fee/K<limit
                margins.append(limit-fee/K)
            boundary=K*limit
            assert abs(K-boundary-(1-beta)*((1+sigma/q)*K-boundary))<1e-14
    report['initial_fee_interiority_checks']=len(margins)
    report['min_initial_fee_fraction_slack']=min(margins)
    (P/'checks/neutrality_analysis.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
