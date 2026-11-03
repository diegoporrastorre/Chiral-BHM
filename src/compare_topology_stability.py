"""Compare zero/frequency-resolved topology with the existing phase sweep."""
from pathlib import Path
import json
import numpy as np
from bdg_topology import winding,finite_svd,zero_frequency_interval,bloch

ROOT=Path(__file__).resolve().parents[1]


def run():
    old=json.loads((ROOT/'data/generated/bdg_phi_sweep.json').read_text())
    results=[]
    for gamma in [.35,.1]:
        rows=[]
        for oldrow in old['rows']:
            phi=np.deg2rad(oldrow['degrees'])
            r=winding(1,phi,.5,.4,gamma)
            fine=winding(1,phi,.5,.4,gamma,nk=16384)
            assert r['winding']==fine['winding']
            interval=zero_frequency_interval(1,.5,.4,gamma)
            expected=None if oldrow['degrees']==0 else int(interval[0]<oldrow['degrees']<interval[1])
            assert r['winding']==expected,(oldrow['degrees'],r,expected)
            rows.append(dict(degrees=oldrow['degrees'],**r,
                             bulk_growth=oldrow['bulk']-gamma,
                             obc40_growth=oldrow['obc']['40']['threshold']-gamma,
                             obc80_growth=oldrow['obc']['80']['threshold']-gamma))
        selected=[]
        for deg in [40,60,70,80,85,90]:
            selected.append(dict(degrees=deg,svd=[finite_svd(N,1,np.deg2rad(deg),.5,.4,gamma)
                                               for N in [20,40,80,120]]))
        # Nonzero frequency example at original phase, omega0=1.5J.
        finite_frequency=dict(omega0=1.5,phi_degrees=60,
            **winding(1,np.pi/3,.5,.4,gamma,omega0=1.5),
            svd=[finite_svd(N,1,np.pi/3,.5,.4,gamma,omega0=1.5) for N in [20,40,80]])
        results.append(dict(gamma=gamma,kappa=2*gamma,interval=interval,rows=rows,
                            selected=selected,finite_frequency=finite_frequency))
        print('LOSS',gamma,'INTERVAL',interval,'FINITE FREQ WIND',finite_frequency['winding'],flush=True)
        for r in selected:
            print(r['degrees'],[(s['N'],s['smin']) for s in r['svd']],flush=True)
    output=dict(parameters=dict(J=1,delta=.5,g=.4,g_c=0),reference_frequency=0,
                conventions='Q=omega0 I-H; k:-pi to pi; forward hopping J exp(+i phi)',results=results)
    (ROOT/'data/generated/bdg_topology_stability.json').write_text(json.dumps(output,indent=2))
    return output


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
    for row,(ax,svax) in zip(result['results'],axes):
        data=row['rows'];x=[r['degrees'] for r in data]
        ax.axvspan(*row['interval'],color='#ddbd55',alpha=.3,label=r'$\nu(0)=+1$')
        ax.plot(x,[r['bulk_growth'] for r in data],'k-',label='Bulk saddle')
        ax.plot(x,[r['obc80_growth'] for r in data],'o-',ms=3,label='OBC N=80')
        ax.axhline(0,color='gray',lw=.8)
        ax.set(title=fr'$\kappa/J={row["kappa"]:g}$; zero-frequency topology',xlabel=r'$\phi$ (degrees)',ylabel='Growth rate / J',xlim=(0,90))
        ax.legend(fontsize=8);ax.grid(alpha=.2)
        for sel in row['selected']:
            if sel['degrees'] not in [60,70,80,85]: continue
            ss=sel['svd'];svax.semilogy([s['N'] for s in ss],[s['smin'] for s in ss],'o-',ms=3,label=fr'$\phi={sel["degrees"]}^\circ$')
        svax.set(title=r'OBC smallest singular value of $-\mathcal{H}$',xlabel='N',ylabel=r'$s_{\min}/J$')
        svax.legend(fontsize=8);svax.grid(alpha=.2)
    fig.suptitle(r'$J=1,\ \tilde\Delta=0.5,\ |g_s|=0.4,\ g_c=0$; shaded regions have nonzero winding')
    fig.savefig(ROOT/'figures/bdg_topology_stability.png',dpi=180)
    plt.close(fig)


if __name__=='__main__':
    plot(run())
