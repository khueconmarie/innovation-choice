"""Exact joint arrival/quality comparisons, admissibility and fee overlap."""
from pathlib import Path
import json,math,csv
import numpy as np
import sympy as sp
from scipy.optimize import brentq
from dynamic_analysis import wtp,H_closed,BETA,BASE,GAIN,LOW,HIGH

P=Path(__file__).resolve().parent

def burden(gamma,T,R,beta=BETA):
    return -math.log1p(-wtp(gamma,R,beta=beta,T=T))

def main():
    a,b=sp.symbols('a b',positive=True);t=sp.symbols('t',integer=True,positive=True)
    h=lambda T:(1-a**T)/(1-a)+a**T/(1-b)
    assert sp.simplify(h(t-1)-h(t)-a**(t-1)*(b-a)/(1-b))==0
    lo=max(wtp(.5,LOW),wtp(2,HIGH))
    hi=min(wtp(.5,HIGH),wtp(2,LOW))
    assert lo<.0839<hi
    cap=BETA**(-2)
    threshold=brentq(lambda R:wtp(.5,R,T=9)-wtp(.5,LOW,T=10),LOW,cap*(1-1e-10))
    comparisons=[]
    for T,Rnew in [(5,HIGH),(10,HIGH),(10,threshold),(10,1.07)]:
        old=wtp(.5,LOW,T=T);same=wtp(.5,LOW,T=T-1);new=wtp(.5,Rnew,T=T-1)
        change=burden(.5,T-1,Rnew)-burden(.5,T,LOW)
        quality=burden(.5,T-1,Rnew)-burden(.5,T-1,LOW)
        arrival=burden(.5,T-1,LOW)-burden(.5,T,LOW)
        assert abs(change-quality-arrival)<1e-14
        comparisons.append(dict(T_old=T,T_new=T-1,R_old=LOW,R_new=Rnew,p_old=old,
                                p_arrival_only=same,p_new=new,quality_component=quality,
                                arrival_component=arrival,total_log_burden_change=change))
    assert max(r['p_new'] for r in comparisons if r['T_old']==5)<wtp(.5,LOW,T=5)
    p4sup=1-GAIN**(-4)
    assert p4sup<wtp(.5,LOW,T=5)
    p9sup=1-GAIN**(-9)
    assert p9sup>wtp(.5,LOW,T=10)
    # The arrival component is a change in both opportunity paths, not a
    # universally negative mechanical truncation of the innovation's benefits.
    earlier_counterexample=dict(gamma=.5,R=1.06,T_old=100,T_new=99,
                               p_old=wtp(.5,1.06,T=100),p_new=wtp(.5,1.06,T=99))
    assert earlier_counterexample['p_new']>earlier_counterexample['p_old']
    rng=np.random.default_rng(380916);maxerror=0.
    for i in range(180):
        gamma=float(rng.uniform(.3,2.7));beta=float(rng.uniform(.85,.97))
        if abs(gamma-1)<.05:gamma=1.1
        T=int(rng.integers(2,40));R0=float(rng.uniform(BASE*GAIN,1.06));R1=R0+.005
        if beta*R1**(1-gamma)>=1:continue
        h0=H_closed(gamma,BASE,R0,beta,T)
        h1=H_closed(gamma,BASE,R0,beta,T-1)
        nu=(1-gamma)/gamma;aa=beta**(1/gamma)*BASE**nu;bb=beta**(1/gamma)*R0**nu
        maxerror=max(maxerror,abs(h1-h0-aa**(T-1)*(bb-aa)/(1-bb)))
        d=burden(gamma,T-1,R1,beta)-burden(gamma,T,R0,beta)
        parts=burden(gamma,T-1,R1,beta)-burden(gamma,T-1,R0,beta)+burden(gamma,T-1,R0,beta)-burden(gamma,T,R0,beta)
        assert abs(d-parts)<1e-13
    assert maxerror<2e-12
    inadmissible=False
    try:wtp(.5,HIGH,beta=.98)
    except (AssertionError,ValueError):inadmissible=True
    assert inadmissible
    report=dict(fee_overlap=[lo,hi],fee_overlap_width_pp=100*(hi-lo),
                beta_ceiling_gamma_half_R106=HIGH**(-.5),R_ceiling_beta096_gamma_half=cap,
                invalid_beta098_rejected=inadmissible,symbolic_H_arrival_difference='pass',
                max_H_difference_error=maxerror,comparisons=comparisons,
                T4_wtp_supremum=p4sup,T9_wtp_supremum=p9sup,compensating_R_T10_to9=threshold,
                earlier_arrival_need_not_reduce_wtp=earlier_counterexample)
    (P/'checks/joint_analysis.json').write_text(json.dumps(report,indent=2)+'\n')
    with (P/'generated/joint_examples.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=comparisons[0].keys());w.writeheader();w.writerows(comparisons)
    lines=[r'\begin{tabular}{rrrrr}',r'\toprule',
           r'$T\to T-1$ & New $R$ & Original $p^*$ & New $p^*$ & Change (pp)\\',r'\midrule']
    for r in comparisons:
        lines.append(f"{r['T_old']} $\\to$ {r['T_new']} & {r['R_new']:.6f} & {100*r['p_old']:.4f}\\% & {100*r['p_new']:.4f}\\% & {100*(r['p_new']-r['p_old']):+.4f}"+r'\\')
    lines.extend([r'\bottomrule',r'\end{tabular}',''])
    (P/'generated/joint_table.tex').write_text('\n'.join(lines))
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
