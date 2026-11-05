"""Add bulk zero-frequency point-gap winding to the saved N=50 stability map."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from compensated_stability_map import bounds, drive, growth as get_growth, CUSP
from scipy.optimize import brentq


def main(plot_only=False):
    folder = Path('data/generated/compensation_n50')
    folder.mkdir(parents=True, exist_ok=True)
    if not (folder/'compensated_stability_map.npz').exists():
        edges=np.linspace(0,CUSP,321); phases=(edges[:-1]+edges[1:])/2
        base=np.sin(np.linspace(0,np.pi/2,81))**2
        extra=np.array([1e-7,1e-6,1e-5,1e-4,.001,.003])
        frac=np.unique(np.r_[base,extra,1-extra]);frac=frac[(frac>0)&(frac<1)]
        low,high=bounds(phases); dens=low[:,None]+(high-low)[:,None]*frac
        rates=np.empty_like(dens); curves=[]
        for i,p in enumerate(phases):
            rates[i]=[get_growth(float(v),float(p)) for v in dens[i]]
            for j in range(len(frac)-1):
                if rates[i,j]*rates[i,j+1]<0:
                    root=brentq(lambda v:get_growth(v,p),dens[i,j],dens[i,j+1],xtol=1e-7)
                    curves.append(dict(phi=float(p),density=float(root),epsilon=float(drive(root,p))))
            if i%20==0: print('Regenerated row',i,'/320',flush=True)
        np.savez_compressed(folder/'compensated_stability_map.npz',phi=phases,density=dens,epsilon=drive(dens,phases[:,None]),growth=rates)
        (folder/'compensated_stability_map.json').write_text(json.dumps(dict(boundary_points=curves,provenance='Regenerated grid; original precision checks are documented in the existing report, not rerun here.')))
    src = np.load(folder/'compensated_stability_map.npz')
    phi = np.broadcast_to(src['phi'][:, None], src['density'].shape)
    n, eps, growth = src['density'], src['epsilon'], src['growth']
    d, g, gamma = .4998 - .0008*n, .0004*n, .5
    x0 = g*g-gamma*gamma-(d+2*np.cos(phi))**2
    xp = g*g-gamma*gamma-(d-2*np.cos(phi))**2
    # Imag det Q = -4 gamma sin(phi) sin(k). In this positive-phi
    # existence lobe the two real-axis crossings determine the winding.
    nu = ((np.sign(xp)-np.sign(x0))/2).astype(int)
    assert np.all(x0 > 0)  # intermediate cubic branch has negative slope
    checks = json.loads((folder/'compensated_topology_map.json').read_text())['numerical_checks'] if plot_only else []
    for nk in (() if plot_only else (4096, 8192)):
        k = np.linspace(-np.pi, np.pi, nk, endpoint=False)
        numerical = np.empty(n.size, int)
        maxstep = np.empty(n.size)
        for start in range(0, n.size, 64):
            sl = slice(start, start+64)
            p, dd, gg = phi.ravel()[sl,None], d.ravel()[sl,None], g.ravel()[sl,None]
            a = 2*np.sin(p)*np.sin(k)
            det = (a-1j*gamma)**2-(dd+2*np.cos(p)*np.cos(k))**2+gg**2
            steps = np.angle(np.roll(det,-1,axis=1)/det)
            numerical[sl] = np.rint(steps.sum(axis=1)/(2*np.pi)).astype(int)
            maxstep[sl] = np.abs(steps).max(axis=1)
        assert np.array_equal(numerical, nu.ravel())
        checks.append(dict(nk=nk, mismatches=0, max_phase_step=float(maxstep.max())))
        print(checks[-1], flush=True)
    np.savez_compressed(folder/'compensated_topology_map.npz', phi=phi, epsilon=eps,
                        density=n, growth=growth, winding=nu, X0=x0, Xpi=xp)
    summary = dict(reference_frequency=0, samples=int(n.size), numerical_checks=checks,
                   min_X0=float(x0.min()), max_Xpi=float(xp.max()),
                   counts={f'{label}_nu_{w}': int(((growth < -1e-7 if label == 'stable' else growth > 1e-7)&(nu==w)).sum())
                           for label in ['stable','unstable'] for w in [-1,0]})
    (folder/'compensated_topology_map.json').write_text(json.dumps(summary,indent=2))
    meta = json.loads((folder/'compensated_stability_map.json').read_text())
    curves = meta['boundary_points']
    tx=np.r_[eps.ravel(),[r['epsilon'] for r in curves]]
    ty=np.r_[phi.ravel(),[r['phi'] for r in curves]]
    rate=np.r_[growth.ravel(),np.zeros(len(curves))]
    tri=Triangulation(tx,ty)
    nt=np.r_[n.ravel(),[r['density'] for r in curves]]
    xt=(.0004*nt)**2-.25-(.4998-.0008*nt-2*np.cos(ty))**2
    dp=np.linspace(0,CUSP,2001);lo,hi=bounds(dp)
    lower,upper=drive(hi,dp),drive(lo,dp)
    for name,limits in [('zoom',(30,40,0,np.deg2rad(72))),('full',(0,85,0,np.pi))]:
        fig,axs=plt.subplots(1,2,figsize=(14,6),sharex=True,sharey=True)
        fig.subplots_adjust(left=.06,right=.985,bottom=.12,top=.84,wspace=.08)
        axs[0].tricontourf(tri,rate,levels=[-100,-1e-7,1e-7,100],colors=['#26866d','#edc64b','#cf6b4b'])
        axs[1].tricontourf(tri,xt,levels=[-1000,0,1000],colors=['#6947ac','#c9dfed'])
        for ax in axs:
            ax.set_facecolor('#e7e7e7')
            if xt.min()<0<xt.max():
                ax.tricontour(tri,xt,levels=[0],colors='#39215f',linewidths=1.5,linestyles='--')
            ax.fill_betweenx(dp,0,lower,color='#e7e7e7',zorder=3)
            ax.fill_betweenx(dp,upper,85,color='#e7e7e7',zorder=3)
            ax.plot(lower,dp,'k',lw=.8,zorder=4);ax.plot(upper,dp,'k',lw=.8,zorder=4)
            ax.axvline(38.51,color='#555',ls=':',lw=.8)
            ax.set(xlim=limits[:2],ylim=limits[2:],xlabel=r'Drive $\epsilon/J$')
            ax.set_yticks(np.arange(7)*np.pi/6,['0',r'$\pi/6$',r'$\pi/3$',r'$\pi/2$',r'$2\pi/3$',r'$5\pi/6$',r'$\pi$'])
            ax.set_ylim(limits[2:])
            ax.grid(alpha=.12)
        axs[1].scatter([r['epsilon'] for r in curves],[r['phi'] for r in curves],s=1,color='#222',zorder=4)
        axs[1].plot([lower[0],upper[0]],[0,0],color='#d73583',lw=3,zorder=5)
        axs[0].set_title('OBC dynamical stability',fontsize=12)
        axs[1].set_title(r'Bulk winding at $\omega_0=0$; black dots: stability boundary',fontsize=11)
        axs[0].set_ylabel(r'Hopping phase $\phi$')
        axs[0].legend(handles=[Patch(color='#26866d',label='Stable'),Patch(color='#cf6b4b',label='Unstable'),Patch(color='#e7e7e7',label='No intermediate root')],loc='upper right')
        axs[1].legend(handles=[Patch(color='#6947ac',label=r'Topological: $\nu=-1$'),Line2D([0],[0],color='#d73583',lw=3,label=r'$\phi=0$: gapless, $\nu$ undefined'),Line2D([0],[0],color='k',lw=.8,label='Root-merger curves: gapless')],loc='upper right',fontsize=9)
        if np.all(nu==-1):
            axs[1].text(.04,.06,'All sampled intermediate states: '+r'$\nu=-1$'+'\nBoth stable and unstable OBC regions are topological.',transform=axs[1].transAxes,fontsize=10,bbox=dict(facecolor='white',alpha=.95,edgecolor='none'),zorder=6)
        fig.suptitle(r'Compensated homogeneous intermediate state | $N=50$, $V=0$, $\Delta=0.5$, $J=\kappa=1$',fontsize=13,y=.96)
        fig.savefig(f'figures/compensated_stability_topology_{name}.png',dpi=180)
        fig.savefig(f'figures/compensated_stability_topology_{name}.svg')
        plt.close(fig)
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    main()
