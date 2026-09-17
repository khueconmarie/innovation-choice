"""Exact dynamic complementarity, physical-path checks, and figure generation.

CRRA utility and deterministic linear accumulation. All numerical inputs are
chosen examples. No empirical calibration or numerical optimizer is used.
"""
from pathlib import Path
import math, json, csv
import numpy as np
import sympy as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P=Path(__file__).resolve().parent
BETA=.96
BASE=1.02
GAIN=1.02
T=5
LOW=1.041
HIGH=1.06
FEE=.0839

def geometric_sum(a,T):
    return sum(a**t for t in range(T))

def H_closed(gamma,ret,common,beta=BETA,T=T):
    nu=(1-gamma)/gamma
    a=beta**(1/gamma)*ret**nu
    b=beta**(1/gamma)*common**nu
    if not (0<beta<1 and gamma>0 and b<1):
        raise ValueError('The infinite-horizon CRRA value-existence condition fails.')
    return geometric_sum(a,T)+a**T/(1-b)

def wtp(gamma,common,beta=BETA,T=T):
    if abs(gamma-1)<1e-10:
        return -math.expm1(-math.log(GAIN)*sum(beta**(s+1) for s in range(T)))
    hb=H_closed(gamma,BASE,common,beta,T)
    hi=H_closed(gamma,BASE*GAIN,common,beta,T)
    return -math.expm1(gamma/(1-gamma)*math.log(hb/hi))

def physical_value(gamma,common,project=False,fee=0.,N=6000):
    # Euler consumption ratios plus the PV budget determine a global optimum.
    # Then capital is reconstructed from the remaining future expenditure to
    # avoid unstable subtraction of nearly exhausted terminal wealth.
    t=np.arange(N,dtype=float)
    rr=np.where(t<T, BASE*(GAIN if project else 1.), common)
    logs=np.r_[0., np.cumsum(np.log(rr[:-1]))]
    nu=(1-gamma)/gamma
    terms=np.exp(t*math.log(BETA)/gamma+nu*logs)
    H=H_closed(gamma,BASE*(GAIN if project else 1.),common)
    c=(1-fee)/H*np.exp((t*math.log(BETA)+logs)/gamma)
    # Only dates 0:200 enter relative feasibility checks; summation is much longer.
    tail=np.cumsum(terms[::-1])[::-1]
    k=np.exp(logs)*(1-fee)/H*tail
    inv=k-c
    # The PV formula gives post-fee resources at date zero, not installed K0.
    k[0]=1.
    initial_budget_error=abs(c[0]+inv[0]+fee-k[0])
    assert initial_budget_error<1e-12
    physical_error=float(np.max(np.abs(k[1:201]-rr[:200]*inv[:200])/
                                  np.maximum(k[1:201],1e-200)))
    physical_error=max(physical_error,initial_budget_error)
    assert np.all(inv[:201]>0)
    euler=float(np.max(np.abs((c[1:201]/c[:200])**gamma-BETA*rr[:200])))
    if abs(gamma-1)<1e-10:
        value=float((BETA**t)@np.log(c))
    else:
        value=float((BETA**t)@(c**(1-gamma)/(1-gamma)))
    pv=float(np.exp(-logs)@c)
    return value,physical_error,euler,abs(pv-(1-fee))

def main():
    for d in ['checks','generated','figures']:(P/d).mkdir(exist_ok=True)
    # Symbolic derivative of exact willingness to pay after a shared later gain.
    gam, z, hb, hi, tail = sp.symbols('gam z hb hi tail',positive=True)
    nu=(1-gam)/gam
    p=1-(hb/hi)**(gam/(1-gam))
    chain=sp.diff(p,hb)*nu*tail+sp.diff(p,hi)*nu*z*tail
    target=(1-p)*tail*(z*hb-hi)/(hb*hi)
    assert sp.simplify(chain-target)==0

    rows=[]; errors=[]
    for gamma in [.5,1.,2.]:
        for common in [LOW,HIGH]:
            pstar=wtp(gamma,common)
            vb,phys0,eul0,bud0=physical_value(gamma,common)
            vi,phys1,eul1,bud1=physical_value(gamma,common,True,FEE)
            vc,phys2,eul2,bud2=physical_value(gamma,common,True,pstar)
            errors.append(dict(gamma=gamma,common=common,indifference=abs(vc-vb),
                               physical=max(phys0,phys1,phys2),euler=max(eul0,eul1,eul2),
                               budget=max(bud0,bud1,bud2)))
            assert (vi>vb)==(pstar>FEE)
            rows.append(dict(gamma=gamma,EIS=1/gamma,common_return=common,
                             willingness_to_pay_share=pstar,fee_share=FEE,
                             net_value=vi-vb,fund=pstar>FEE,
                             optimal_initial_investment_share=1-1/H_closed(gamma,BASE*GAIN,common),
                             limiting_consumption_growth=math.log(BETA*common)/gamma,
                             innovation_growth_attribution=0.))
    assert [r['fund'] for r in rows]==[False,True,True,True,True,False]
    assert max(e['indifference'] for e in errors)<1e-10
    assert max(e['physical'] for e in errors)<1e-11

    # General finite-horizon paths; arbitrary strictly positive returns and
    # nonnegative early project gains, not just constant-return examples.
    rng=np.random.default_rng(370916)
    derivative_error=0.; sign_checks=0; exposure_error=0.
    for j in range(240):
        gamma=float(rng.uniform(.3,2.8))
        if abs(gamma-1)<.05:gamma+=.1
        beta=float(rng.uniform(.86,.98)); n=120
        T0=int(rng.integers(2,20)); k=int(rng.integers(T0,55))
        rb=np.exp(rng.uniform(-.02,.05,n-1))
        ri=rb.copy();ri[:T0]*=np.exp(rng.uniform(.001,.04,T0))
        def finite_H(r):
            lg=np.r_[0.,np.cumsum(np.log(r))]
            return np.exp(np.arange(n)*math.log(beta)/gamma+(1-gamma)/gamma*lg)
        wb=finite_H(rb);wi=finite_H(ri)
        HB,HI=wb.sum(),wi.sum();Z=wi[-1]/wb[-1]
        p0=1-(HB/HI)**(gamma/(1-gamma))
        deriv=(1-p0)*wb[k+1:].sum()*(Z*HB-HI)/(HB*HI)
        step=2e-5
        pp=[]
        for eps in [-step,step]:
            rbe=rb.copy();rie=ri.copy();rbe[k]*=math.exp(eps);rie[k]*=math.exp(eps)
            pp.append(1-(finite_H(rbe).sum()/finite_H(rie).sum())**(gamma/(1-gamma)))
        err=abs((pp[1]-pp[0])/(2*step)-deriv)
        derivative_error=max(derivative_error,err)
        assert math.copysign(1,deriv)==math.copysign(1,1-gamma)
        sign_checks+=1
        # Independent finite-horizon value derivative equals beta^k u'(c_k) I_k.
        r=ri; ww=wi; H=HI; lg=np.r_[0.,np.cumsum(np.log(r))]
        c=np.exp((np.arange(n)*math.log(beta)+lg)/gamma)/H
        capital=np.exp(lg)*np.cumsum(ww[::-1])[::-1]/H
        investment=capital-c
        weight=beta**k*c[k]**(-gamma)*investment[k]
        indirect=H**(gamma-1)*ww[k+1:].sum()
        exposure_error=max(exposure_error,abs(weight-indirect))
        assert np.all(investment[:-1]>0)
    assert derivative_error<1e-8
    assert exposure_error<1e-10

    lines=[r'\begin{tabular}{rrrrll}',r'\toprule',
           r'$\gamma$ & EIS & Low substitute & High substitute & Low: fund? & High: fund?\\',r'\midrule']
    for i in range(0,len(rows),2):
        l,h=rows[i:i+2]
        lines.append(f"{l['gamma']:.1f} & {l['EIS']:.1f} & {100*l['willingness_to_pay_share']:.4f}\\% & "
                     f"{100*h['willingness_to_pay_share']:.4f}\\% & "
                     f"{'Yes' if l['fund'] else 'No'} & {'Yes' if h['fund'] else 'No'} "+r'\\')
    lines.extend([r'\bottomrule',r'\end{tabular}',''])
    (P/'generated/dynamic_table.tex').write_text('\n'.join(lines))
    with (P/'generated/dynamic_examples.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)

    # The publication figure is generated separately by render_figure.py.
    

    report={'model':'Deterministic linear accumulation, CRRA, endogenous optimal consumption/investment',
            'inputs':{'beta':BETA,'base_return':BASE,'innovation_multiplier':GAIN,'obsolescence_date':T,
                      'low_common_return':LOW,'high_common_return':HIGH,'project_fee_share':FEE},
            'examples':rows,'symbolic_WTP_derivative':'pass','physical_checks':errors,
            'common_physical_initial_capital_checked':True,
            'random_path_sign_checks':sign_checks,'central_difference_step':2e-5,
            'derivative_max_error':derivative_error,'investment_exposure_max_error':exposure_error,
            'interpretation':'Theorem has a global optimum; physical checks verify its recurrences, not an independent optimizer.'}
    (P/'checks/dynamic_analysis.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'examples':rows,'checks':{'random_paths':sign_checks,'derivative_error':derivative_error,
                                             'exposure_error':exposure_error,'physical':errors}},indent=2))

if __name__=='__main__':main()
