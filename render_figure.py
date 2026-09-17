"""Journal figure with line styles that remain distinguishable in grayscale."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dynamic_analysis import wtp, H_closed, BASE, GAIN, LOW, HIGH, FEE

P=Path(__file__).resolve().parent
(P/'figures').mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'serif','font.size':10,'pdf.fonttype':42,
    'axes.spines.top':False,'axes.spines.right':False})
fig,axs=plt.subplots(1,2,figsize=(10,3.8))
xx=np.linspace(LOW,HIGH,160)
for gamma,color,style in [(.5,'#246483','-'),(1.,'#4d4d4d','--'),(2.,'#9d3f30','-.')]:
    axs[0].plot(xx,[100*wtp(gamma,x) for x in xx],color=color,ls=style,
                label=f'EIS = {1/gamma:g}',lw=1.9)
    axs[1].plot(xx,[1-1/H_closed(gamma,BASE*GAIN,x) for x in xx],
                color=color,ls=style,lw=1.9)
axs[0].axhline(100*FEE,color='black',lw=1.,ls=':',label='Project fee')
axs[0].set(xlabel='Common gross return after obsolescence',ylabel='Maximum initial fee (% of resources)')
axs[1].set(xlabel='Common gross return after obsolescence',ylabel='Initial investment / post-fee resources')
for i,ax in enumerate(axs):
    ax.set_xticks([1.041,1.047,1.053,1.060])
    ax.set_xticklabels(['1.041','1.047','1.053','1.060'])
    ax.text(.02,.98,f'({chr(97+i)})',transform=ax.transAxes,va='top')
fig.legend(*axs[0].get_legend_handles_labels(),loc='upper center',ncol=4,
           frameon=False,bbox_to_anchor=(.5,1.02),fontsize=9)
fig.tight_layout(rect=[0,0,1,.91])
fig.savefig(P/'figures/Fig1.pdf',bbox_inches='tight')
fig.savefig(P/'figures/Fig1.png',bbox_inches='tight',dpi=220)
plt.close(fig)
