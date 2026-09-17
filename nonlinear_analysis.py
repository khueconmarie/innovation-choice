"""One-state concave investment problem with feasible and affine value bounds.

Capital is the only endogenous state. A concave prefix objective plus a
stationary affine supersolution supplies the upper bound. Holding terminal
capital fixed supplies an implementable infinite continuation and lower bound.
The residual correction uses a global feasible capital box, not a local Hessian.
"""
from dataclasses import dataclass,replace
from pathlib import Path
import json, math,csv
import numpy as np
from scipy.linalg import solve_banded
from scipy.optimize import brentq

P=Path(__file__).resolve().parent

@dataclass(frozen=True)
class Economy:
    eta:float=.36
    sigma:float=.90
    gamma:float=.5
    beta:float=.96
    qbase:float=1.02
    gain:float=1.02
    T:int=5
    K0:float=1.
    N:int=160

    @property
    def A(self):
        # The common pre-innovation steady state is K0=1, for every curvature.
        return (1/self.beta-self.sigma)/(self.eta*self.qbase)

    def u(self,c):
        return np.log(c) if self.gamma==1 else c**(1-self.gamma)/(1-self.gamma)

    def solve(self,common,access=False,fee=0.,initial=None):
        n=self.N;beta=self.beta;sig=self.sigma;eta=self.eta;g=self.gamma
        q=np.full(n,common);q[:self.T]=self.qbase*(self.gain if access else 1.)
        ks=(common*self.A*eta/(1/beta-sig))**(1/(1-eta))
        cs=self.A*ks**eta-(1-sig)*ks/common
        slope=cs**(-g)/(beta*common)
        vs=float(self.u(cs)/(1-beta))
        kmax=max(self.K0,(max(q.max(),common)*self.A/(1-sig))**(1/(1-eta)))
        discount=beta**np.arange(n)

        def allocation(x):
            k=np.r_[self.K0,x]
            y=self.A*k[:-1]**eta
            inv=(k[1:]-sig*k[:-1])/q
            c=y-inv;c[0]-=fee
            return k,c,inv

        if initial is None:
            share=(1-sig)*ks/(common*self.A*ks**eta)
            k=[self.K0]
            for t in range(n):
                resources=self.A*k[-1]**eta-(fee if t==0 else 0.)
                k.append(sig*k[-1]+q[t]*share*resources)
            x=np.array(k[1:])
        else:x=initial.copy()

        def system(x,jac=True):
            k,c,inv=allocation(x)
            if np.min(c)<=0 or np.min(inv)<=0 or np.min(k)<=0:return None
            fp=self.A*eta*k[:-1]**(eta-1)
            a=fp+sig/q
            res=np.empty(n)
            res[:-1]=-g*np.log(c[:-1])-np.log(q[:-1])-math.log(beta)+g*np.log(c[1:])-np.log(a[1:])
            res[-1]=-g*math.log(c[-1])-math.log(q[-1])-math.log(beta)-math.log(slope)
            if not jac:return res
            fpp=self.A*eta*(eta-1)*k[1:-1]**(eta-2)
            diag=g/(q*c)
            diag[:-1]+=g*a[1:]/c[1:]-fpp/a[1:]
            upper=-g/(q[1:]*c[1:])
            lower=-g*a[1:]/c[1:]
            band=np.zeros((3,n));band[1]=diag;band[0,1:]=upper;band[2,:-1]=lower
            return res,band

        for iteration in range(100):
            result=system(x)
            if result is None:raise RuntimeError('infeasible initial path')
            res,band=result;norm=float(np.max(np.abs(res)))
            if norm<2e-12:break
            step=solve_banded((1,1),band,-res)
            scale=1.
            for _ in range(60):
                trial=x+scale*step
                if np.all(trial>0) and np.all(trial<kmax):
                    rr=system(trial,False)
                    if rr is not None and np.max(np.abs(rr))<norm*(1-1e-4*scale):break
                scale*=.5
            else:raise RuntimeError(('line search',norm))
            x=trial
        else:raise RuntimeError(('Newton nonconvergence',norm))
        k,c,inv=allocation(x)
        prefix=float(discount@self.u(c))
        affine=vs+slope*(k[-1]-ks)
        holdc=self.A*k[-1]**eta-(1-sig)*k[-1]/common
        assert holdc>0 and np.min(inv)>0
        lower=prefix+beta**n*float(self.u(holdc)/(1-beta))
        upper_at_x=prefix+beta**n*affine
        # Full gradient of the concave prefix objective, in initial utility units.
        fp=self.A*eta*k[:-1]**(eta-1)
        grad=-discount*c**(-g)/q
        grad[:-1]+=discount[1:]*c[1:]**(-g)*(fp[1:]+sig/q[1:])
        grad[-1]+=beta**n*slope
        correction=float(np.maximum(grad*(kmax-x),-grad*x).sum())
        # A stated numerical rounding allowance, separate from the analytic bound.
        pad=2e-11*(1+abs(lower)+abs(upper_at_x))
        upper=upper_at_x+correction+pad;lower-=pad
        assert lower<=upper
        physical=float(np.max(abs(k[1:]-sig*k[:-1]-q*inv)))
        budget=float(np.max(abs(c+inv+np.r_[fee,np.zeros(n-1)]-self.A*k[:-1]**eta)))
        return dict(lower=lower,upper=upper,mid=(lower+upper)/2,
                    bound_width=upper-lower,gradient_box_correction=correction,
                    rounding_allowance=pad,euler_log_residual=norm,
                    physical_residual=physical,budget_residual=budget,
                    minimum_investment=float(inv.min()),initial_investment_share=float(inv[0]/(self.A-fee)),
                    terminal_capital=float(k[-1]),steady_capital=ks,capital_box=kmax,
                    iterations=iteration,path=k,consumption=c,investment=inv)

    def wtp(self,common):
        outside=self.solve(common)
        evaluated=[]
        def difference(f):
            s=self.solve(common,True,f*self.A)
            evaluated.append((f,s))
            return s['mid']-outside['mid']
        root=brentq(difference,0.,.15,xtol=2e-13,rtol=1e-13)
        at=self.solve(common,True,root*self.A)
        # Two endpoint decisions bracket the infinite-horizon compensating fee.
        marginal=self.A*at['consumption'][0]**(-self.gamma)
        error=(at['bound_width']+outside['bound_width'])/marginal
        half=max(2e-9,2*error)
        lo=self.solve(common,True,(root-half)*self.A)
        hi=self.solve(common,True,(root+half)*self.A)
        assert lo['lower']>outside['upper'] and hi['upper']<outside['lower']
        summary={k:v for k,v in at.items() if k not in ('path','consumption','investment')}
        return dict(eta=self.eta,sigma=self.sigma,gamma=self.gamma,EIS=1/self.gamma,N=self.N,
                    A=self.A,common_quality=common,fee_share_output=root,
                    fee_interval=[root-half,root+half],outside_value_interval=[outside['lower'],outside['upper']],
                    access_value_interval=[at['lower'],at['upper']],
                    lower_fee_value_margin=lo['lower']-outside['upper'],
                    upper_fee_value_margin=hi['upper']-outside['lower'],checks=summary)

def main():
    cases=[];extended=[];envelopes=[];benchmarks=[];paths={}
    for eta,sigma in [(.90,.20),(.90,.90),(.36,.20),(.36,.90)]:
        # High survival and weak diminishing returns yield a slow transition.
        # Use a longer prefix there to make the value brackets informative.
        n=480 if (eta,sigma)==(.90,.90) else 160
        n_long=640 if n==480 else 240
        for gamma in [.5,1.,2.]:
            economy=Economy(eta=eta,sigma=sigma,gamma=gamma,N=n)
            for common in [1.041,1.06]:
                result=economy.wtp(common);cases.append(result)
                long=replace(economy,N=n_long).wtp(common);extended.append(long)
                assert result['fee_interval'][0]<long['fee_share_output']<result['fee_interval'][1]
                path=economy.solve(common,True,result['fee_share_output']*economy.A)
                key=f'eta{eta}_sigma{sigma}_gamma{gamma}_q{common}'
                for variable in ['path','consumption','investment']:paths[key+'_'+variable]=path[variable]
            # Full quality perturbation after T, holding the compensated fee
            # fixed for the envelope comparison. This is Proposition 1 evaluated
            # using physical investment and consumption, not the WTP formula.
            e=replace(economy,N=n_long);q=1.05;fee=e.wtp(q)['fee_share_output']
            access=e.solve(q,True,fee*e.A);outside=e.solve(q)
            exposure=e.beta**np.arange(e.N)*(access['consumption']**(-gamma)*access['investment']
                       -outside['consumption']**(-gamma)*outside['investment'])
            predicted=float(exposure[e.T:].sum()/(e.A*access['consumption'][0]**(-gamma)))
            step=1e-4
            fd=(e.wtp(q*math.exp(step))['fee_share_output']-e.wtp(q*math.exp(-step))['fee_share_output'])/(2*step)
            envelopes.append(dict(eta=eta,sigma=sigma,gamma=gamma,N=e.N,quality=q,
                                  envelope=predicted,finite_difference=fd,step=step,error=abs(predicted-fd)))
            if n==480:
                # Check sensitivity to step size in the slowly converging cell;
                # the finite-difference error is separate from value brackets.
                perturbations=[]
                for h in [1e-3,5e-4]:
                    other=(e.wtp(q*math.exp(h))['fee_share_output']-e.wtp(q*math.exp(-h))['fee_share_output'])/(2*h)
                    perturbations.append(dict(step=h,finite_difference=other,error=abs(predicted-other)))
                envelopes[-1]['additional_steps']=perturbations
    for eta in [.36,.9]:
        e=Economy(eta=eta,sigma=0.,gamma=1.)
        exact=-math.expm1(-eta*e.beta*sum(e.beta**t for t in range(e.T))*math.log(e.gain))
        for q in [1.041,1.06]:
            computed=e.wtp(q)['fee_share_output']
            benchmarks.append(dict(eta=eta,q=q,exact=exact,computed=computed,error=abs(exact-computed)))
    changes=[]
    lines=[r'\begin{tabular}{rrrrrr}',r'\toprule',
           r'$\eta$ & $\sigma$ & EIS & Low $q^S$: $F^*/Y_0$ & High $q^S$: $F^*/Y_0$ & Change (pp)\\',r'\midrule']
    for i,(low,high) in enumerate(zip(cases[::2],cases[1::2])):
        delta=high['fee_share_output']-low['fee_share_output']
        enclosure=[high['fee_interval'][0]-low['fee_interval'][1],high['fee_interval'][1]-low['fee_interval'][0]]
        assert enclosure[0]*enclosure[1]>0
        changes.append(dict(eta=low['eta'],sigma=low['sigma'],EIS=low['EIS'],change=delta,change_interval=enclosure))
        if i and i%3==0:lines.append(r'\midrule')
        lines.append(f"{low['eta']:.2f} & {low['sigma']:.2f} & {low['EIS']:.1f} & {100*low['fee_share_output']:.5f}\\% & {100*high['fee_share_output']:.5f}\\% & {100*delta:+.5f}"+r'\\')
    lines.extend([r'\bottomrule',r'\end{tabular}',''])
    (P/'generated/nonlinear_table.tex').write_text('\n'.join(lines))
    np.savez_compressed(P/'generated/nonlinear_paths.npz',**paths)
    with (P/'generated/nonlinear_examples.csv').open('w',newline='') as f:
        keys=['eta','sigma','gamma','EIS','N','A','common_quality','fee_share_output']
        writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader()
        writer.writerows({k:row[k] for k in keys} for row in cases)
    report=dict(cases=cases,extended_horizon=extended,sign_intervals=changes,
                physical_envelope_comparisons=envelopes,exact_log_full_depreciation_checks=benchmarks,
                max_horizon_fee_change=max(abs(a['fee_share_output']-b['fee_share_output']) for a,b in zip(cases,extended)),
                max_fee_interval_width=max(c['fee_interval'][1]-c['fee_interval'][0] for c in cases),
                max_envelope_error=max(e['error'] for e in envelopes),
                max_exact_benchmark_error=max(b['error'] for b in benchmarks),
                bound_scope='Analytic affine upper bound and feasible continuation lower bound, numerically evaluated with gradient-box correction and explicit rounding allowance; not directed-rounding interval arithmetic.')
    report['normalization']='A=(1/beta-sigma)/(eta*qbase) keeps pre-innovation steady capital K0=1 in all four cells; A is not held fixed across cells.'
    by_key={(c['eta'],c['sigma'],c['EIS']):c for c in changes}
    def contrast(after,before):
        return dict(change=after['change']-before['change'],
                    interval=[after['change_interval'][0]-before['change_interval'][1],
                              after['change_interval'][1]-before['change_interval'][0]])
    factorial=[]
    for eis in [2.,1.,.5]:
        for sigma in [.2,.9]:
            factorial.append(dict(EIS=eis,comparison='eta .90 to .36',fixed_sigma=sigma,
                                  **contrast(by_key[(.36,sigma,eis)],by_key[(.9,sigma,eis)])))
        for eta in [.9,.36]:
            factorial.append(dict(EIS=eis,comparison='sigma .20 to .90',fixed_eta=eta,
                                  **contrast(by_key[(eta,.9,eis)],by_key[(eta,.2,eis)])))
    report['factorial_contrasts']=factorial
    report['max_euler_log_residual']=max(c['checks']['euler_log_residual'] for c in cases+extended)
    report['max_physical_residual']=max(c['checks']['physical_residual'] for c in cases+extended)
    report['max_budget_residual']=max(c['checks']['budget_residual'] for c in cases+extended)
    assert report['max_horizon_fee_change']<3e-10
    # The new high-survival cell has larger value-difference noise than the
    # inherited cells; report the observed error, with its step-size check.
    assert report['max_envelope_error']<5e-8
    assert report['max_exact_benchmark_error']<3e-13
    (P/'checks/nonlinear_analysis.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['cases','extended_horizon']},indent=2))

if __name__=='__main__':main()
