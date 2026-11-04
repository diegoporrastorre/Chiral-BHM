"""Plot saved stationary density profiles; no new stationary solves."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/generated/homogeneous_obc_n100'

def main():
    records=json.loads((DATA/'scan.json').read_text())[:161]
    refined=json.loads((DATA/'refinement.json').read_text())
    lookup={round(r['epsilon'],6):r for r in records}
    lookup.update({round(r['epsilon'],6):r for r in refined})
    x=np.arange(1,101)
    groups=[('interval', [35.25,37,38,38.5,39,40,42,44.75]),
            ('crossover',[38.25,38.40,38.50,38.55,38.65,38.75])]
    for name,values in groups:
        fig,axes=plt.subplots(len(values)//2,2,figsize=(12,2.5*len(values)//2),
                              sharex=True,sharey=True,layout='constrained')
        for ax,eps in zip(axes.flat,values):
            r=lookup[eps];n=np.array(r['density']);mask=np.array(r['nu'])!=0
            ax.plot(x,n,color='#265d91',lw=1.6,label='stationary density')
            ids=np.where(mask)[0]
            for g in np.split(ids,np.where(np.diff(ids)>1)[0]+1):
                if len(g):ax.axvspan(x[g[0]]-.5,x[g[-1]]+.5,color='#eaaa34',alpha=.25,lw=0)
            ax.plot(x[mask],n[mask],'.',color='#b05b00',ms=3,label='nonzero local winding')
            ax.set(title=rf'$\epsilon/J={eps:.2f}$',xlim=(.5,100.5),ylim=(0,6500))
            ax.grid(alpha=.15)
        for ax in axes[-1]:ax.set_xlabel('Site j')
        for ax in axes[:,0]:ax.set_ylabel(r'Density $|\alpha_j|^2$')
        handles,labels=axes.flat[0].get_legend_handles_labels()
        fig.legend(handles,labels,loc='outside upper center',ncol=2,frameon=False)
        fig.savefig(ROOT/f'figures/topological_segment_density_{name}.png',dpi=180)
        plt.close(fig)

if __name__=='__main__':main()
