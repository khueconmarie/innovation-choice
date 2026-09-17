"""The valuation comparison and its compensated consumption ordering."""
from pathlib import Path
import csv, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dynamic_analysis import wtp,H_closed,BASE,GAIN,LOW,HIGH,FEE,BETA,T
P=Path(__file__).resolve().parent
for d in ['figures','generated','checks']:(P/d).mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'serif','font.size':10,'pdf.fonttype':42,
    'axes.spines.top':False,'axes.spines.right':False})
fig,axs=plt.subplots(1,2,figsize=(10,3.9))
xx=np.linspace(LOW,HIGH,160);dates=np.arange(16)
curves=[];consumption=[]
for gamma,color,style in [(.5,'#246483','-'),(1.,'#4d4d4d','--'),(2.,'#9d3f30','-.')]:
    payments=np.array([wtp(gamma,x) for x in xx])
    axs[0].plot(xx,100*payments,color=color,ls=style,label=f'EIS = {1/gamma:g}',lw=1.9)
    curves.extend(dict(gamma=gamma,common_return=x,fee_fraction=p) for x,p in zip(xx,payments))
    p=wtp(gamma,HIGH);H0=H_closed(gamma,BASE,HIGH);HI=H_closed(gamma,BASE*GAIN,HIGH)
    log_ratio=np.log1p(-p)+np.log(H0/HI)+np.minimum(dates,T)*np.log(GAIN)/gamma
    assert log_ratio[0]<0 and np.all(log_ratio[dates>=T]>0)
    # Directly reconstruct consumption at the compensated resource amounts.
    g0=BASE**np.minimum(dates,T)*HIGH**np.maximum(dates-T,0)
    gi=g0*GAIN**np.minimum(dates,T)
    cb=(BETA**dates*g0)**(1/gamma)/H0
    ci=(1-p)*(BETA**dates*gi)**(1/gamma)/HI
    assert np.max(np.abs(np.log(ci/cb)-log_ratio))<2e-14
    axs[1].plot(dates,log_ratio,color=color,ls=style,lw=1.9)
    consumption.extend(dict(gamma=gamma,date=int(t),common_return=HIGH,fee_fraction=p,
        C_outside=b,C_access=i,log_ratio=l) for t,b,i,l in zip(dates,cb,ci,log_ratio))
axs[0].axhline(100*FEE,color='black',lw=1.,ls=':',label='Fixed project fee')
axs[0].set(xlabel='Common gross return after date 5',ylabel='Reservation payment (% of resources)')
axs[0].set_xticks([1.041,1.047,1.053,1.060]);axs[0].set_xticklabels(['1.041','1.047','1.053','1.060'])
axs[1].axhline(0,color='black',lw=.8);axs[1].axvline(T,color='#888888',ls=':',lw=1.)
axs[1].set(xlabel='Date',ylabel='Compensated log consumption ratio',xticks=[0,5,10,15])
for i,ax in enumerate(axs):ax.text(.02,.98,f'({chr(97+i)})',transform=ax.transAxes,va='top')
fig.legend(*axs[0].get_legend_handles_labels(),loc='upper center',ncol=4,
           frameon=False,bbox_to_anchor=(.5,1.02),fontsize=9)
fig.tight_layout(rect=[0,0,1,.91])
fig.savefig(P/'figures/Fig1.pdf',bbox_inches='tight');fig.savefig(P/'figures/Fig1.png',bbox_inches='tight',dpi=220)
plt.close(fig)
for name,rows in [('reservation_curves',curves),('compensated_consumption_paths',consumption)]:
    with (P/'generated'/f'{name}.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
(P/'checks/figure_analysis.json').write_text(json.dumps(dict(passed=True,curve_rows=len(curves),
    consumption_rows=len(consumption),right_panel_common_return=HIGH,
    right_panel_fee='Own reservation fee in each preference specification',ordering_checked=True),indent=2)+'\n')
