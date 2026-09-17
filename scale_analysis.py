"""Resource ceilings and finite changes in the linear accumulation economy.

Chosen examples, with independent present-value feasibility checks. The
structural proof is in proofs_linear.tex; these tests do not replace it.
"""
from pathlib import Path
import csv
import json
import math
import numpy as np
from dynamic_analysis import wtp, BETA, BASE, GAIN, LOW, HIGH

P = Path(__file__).resolve().parent


def main():
    for name in ('generated', 'checks'):
        (P/name).mkdir(exist_ok=True)
    rows = []
    for T in [5, 10, 20, 40]:
        low, high = wtp(.5, LOW, T=T), wtp(.5, HIGH, T=T)
        cap = -math.expm1(-T*math.log(GAIN))
        assert 0 < low < high < cap
        rows.append(dict(T=T, payment_low=low, payment_high=high,
                         change=high-low, resource_ceiling=cap,
                         remaining_room=cap-low,
                         relative_change=(high-low)/low))
    with (P/'generated/scale_examples.csv').open('w', newline='') as f:
        writer=csv.DictWriter(f, fieldnames=rows[0]); writer.writeheader(); writer.writerows(rows)
    lines=[r'\begin{tabular}{rrrrr}', r'\toprule',
           r'$T$ & Low payment & High payment & Change (pp) & Remaining room (pp)\\',
           r'\midrule']
    for r in rows:
        lines.append(f"{r['T']} & {100*r['payment_low']:.4f}\\% & "
                     f"{100*r['payment_high']:.4f}\\% & {100*r['change']:.4f} & "
                     f"{100*r['remaining_room']:.4f} "+r'\\')
    lines += [r'\bottomrule', r'\end{tabular}', '']
    (P/'generated/scale_table.tex').write_text('\n'.join(lines))

    rng=np.random.default_rng(390916)
    errors=[]; slacks=[]; distances=[]
    for gamma in [.35, .5, 1., 1.7, 3.]:
        for _ in range(80):
            n=180; T=int(rng.integers(1, 40)); beta=float(rng.uniform(.85,.985))
            r0=np.exp(rng.uniform(-.04,.08,n-1)); ri=r0.copy()
            ri[:T] *= np.exp(rng.uniform(.001,.035,T))
            lg0=np.r_[0.,np.cumsum(np.log(r0))]
            lgi=np.r_[0.,np.cumsum(np.log(ri))]
            L=math.exp(float(np.sum(np.log(ri[:T]/r0[:T]))))
            dates=np.arange(n)
            w0=np.exp(dates*math.log(beta)/gamma+(1-gamma)/gamma*lg0)
            wi=np.exp(dates*math.log(beta)/gamma+(1-gamma)/gamma*lgi)
            H0,HI=float(w0.sum()),float(wi.sum())
            c0=np.exp((dates*math.log(beta)+lg0)/gamma)/H0
            ci=np.exp((dates*math.log(beta)+lgi)/gamma)/(L*HI)
            pvI=float(np.exp(-lgi)@ci)
            pv0=float(np.exp(-lg0)@ci)
            slack=1-pv0
            assert slack>0 and abs(pvI-1/L)<1e-12
            # Map the access allocation into a feasible outside allocation,
            # consuming the positive PV slack at date zero.
            mapped=ci.copy(); mapped[0]+=slack
            residual_pv=np.cumsum((mapped*np.exp(-lg0))[::-1])[::-1]
            capital=np.exp(lg0)*residual_pv
            inv=capital-mapped
            assert np.min(inv[:-1])>=-1e-10
            errors.append(float(np.max(np.abs(capital[1:]-r0*inv[:-1])/
                                       np.maximum(capital[1:],1e-100))))
            if gamma==1:
                value0=float((beta**dates)@np.log(c0))
                valueI=float((beta**dates)@np.log(ci))
                payment=-math.expm1(-float((beta**dates)@(lgi-lg0))/H0)
            else:
                value0=float((beta**dates)@(c0**(1-gamma)/(1-gamma)))
                valueI=float((beta**dates)@(ci**(1-gamma)/(1-gamma)))
                payment=-math.expm1(gamma/(1-gamma)*math.log(H0/HI))
            assert value0>valueI and 0<payment<1-1/L
            slacks.append(slack); distances.append(1-1/L-payment)

    Rmax=BETA**(-2)
    approaches=[]
    for T in [5,20]:
        cap=1-GAIN**(-T); previous=0.
        for epsilon in [1e-3,1e-5,1e-7,1e-9]:
            R=Rmax*(1-epsilon); p=wtp(.5,R,T=T)
            assert previous<p<cap
            approaches.append(dict(T=T,epsilon=epsilon,R=R,payment=p,gap_to_ceiling=cap-p))
            previous=p
    report={'inputs':{'beta':BETA,'gamma':.5,'r0':BASE,'rI':BASE*GAIN,
                      'low_R':LOW,'high_R':HIGH,'durations':[5,10,20,40]},
            'rows':rows,'finite_value_margin_1_minus_b':1-BETA**2*HIGH,
            'random_resource_comparisons':len(slacks),'min_PV_slack':min(slacks),
            'min_payment_ceiling_gap':min(distances),'max_physical_error':max(errors),
            'boundary_approaches':approaches,
            'scope':'Linear accumulation resource ceiling; finite values only; no empirical claim.'}
    assert max(errors)<1e-11
    (P/'checks/scale_analysis.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
