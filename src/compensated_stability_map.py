"""2D compensated-intermediate OBC stability map, N=50, no dynamics."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json
import numpy as np
from scipy.linalg import eigvals
from scipy.optimize import brentq
from bulk_stability import matrix
from compensated_phi_stability import precise
from mean_field import Chain
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib.patches import Patch
from matplotlib.tri import Triangulation

DEST=Path('data/generated/compensation_n50')
N=50;U=-.0002;gamma=.5;u=2*U
FIGURE_PREFIX='compensated_stability_map'
PLOT_RANGES=[('full',(0,85,0,np.pi)),('zoom',(30,40,0,np.deg2rad(72)))]
PLOT_TITLE=None
CUSP=float(np.arccos((np.sqrt(3)*gamma-.4998)/2))

def bounds(phi):
    c=.4998+2*np.cos(phi)
    rad=np.sqrt(np.maximum(c*c-3*gamma*gamma,0))
    return (-2*c+rad)/(3*u),(-2*c-rad)/(3*u)

def drive(n,phi):
    return np.sqrt(n*((.4998+2*np.cos(phi)+u*n)**2+gamma*gamma))

def growth(n,phi):
    return float(eigvals(matrix(N,1,phi,.4998+4*U*n,abs(2*U)*n,gamma),check_finite=False).imag.max())

def main():
    # Cubic derivative is negative between its two turning densities. Sample
    # there directly to avoid missing thin strips next to existence curves.
    phase_edges=np.linspace(0,CUSP,321)
    phi=(phase_edges[:-1]+phase_edges[1:])/2
    base=np.sin(np.linspace(0,np.pi/2,81))**2
    f=np.unique(np.r_[base,[1e-7,1e-6,1e-5,1e-4,.001,.003],1-np.array([1e-7,1e-6,1e-5,1e-4,.001,.003])])
    f=f[(f>0)&(f<1)];f_edges=np.r_[0,(f[:-1]+f[1:])/2,1]
    nlo,nhi=bounds(phi);density=nlo[:,None]+(nhi-nlo)[:,None]*f
    epsilon=drive(density,phi[:,None]);rates=np.full(density.shape,np.nan)
    path=DEST/'compensated_stability_map.npz'
    if path.exists():
        old=np.load(path)
        if np.array_equal(old['phi'],phi) and np.array_equal(old['fraction'],f):rates=old['growth'].copy()
    for i,p in enumerate(phi):
        if np.isfinite(rates[i]).all():continue
        rates[i]=[growth(float(n),float(p)) for n in density[i]]
        if i%20==0 or i==len(phi)-1:
            np.savez_compressed(path,phi=phi,phase_edges=phase_edges,fraction=f,fraction_edges=f_edges,density=density,epsilon=epsilon,growth=rates)
            print('row',i,'of',len(phi),'phi degrees',np.rad2deg(p),flush=True)
    # Resolve zero-growth curves along each row, rather than trusting colored cells.
    curves=[]
    for i,p in enumerate(phi):
        for j in range(len(f)-1):
            if rates[i,j]*rates[i,j+1]<0:
                n=brentq(lambda n:growth(n,p),density[i,j],density[i,j+1],xtol=1e-7)
                curves.append(dict(phi=float(p),epsilon=float(drive(n,p)),density=float(n)))
    selected=[]
    for degrees in [0,20,40,60,75]:
        p=np.deg2rad(degrees);a,b=bounds(p)
        for fraction in [.005,.5,.995]:
            n=a+(b-a)*fraction;e=float(drive(n,p));g=growth(n,p)
            checks=[precise(degrees,precision,epsilon=e,initial_density=float(n),N=N) for precision in [100,160]]
            selected.append(dict(degrees=degrees,epsilon=e,density=float(n),growth=g,checks=checks))
            assert checks[0]['counts']==checks[1]['counts']
            assert (g<0)==(checks[1]['counts']['0']==0)
            # Check the actual complex compensated model, not just gauge-rotated pairing.
            c=.4998+2*np.cos(p);alpha=-e/(c+u*n-.5j)
            m=Chain(N=N,drive=e,phi=p);m.pump=m.pump.astype(complex)
            m.pump[0]+=np.exp(1j*p)*alpha;m.pump[-1]+=np.exp(-1j*p)*alpha
            residual=float(max(abs(m.force(np.full(N,alpha)))))
            assert residual<1e-8;selected[-1]['stationary_residual']=residual
            print('precision',degrees,e,g,checks[-1]['counts']['0'],flush=True)
    summary=dict(N=N,parameters=dict(J=1,U=U,V=0,delta=.5,kappa=1),cusp_phi=CUSP,cusp_degrees=float(np.rad2deg(CUSP)),
                 samples=int(rates.size),phase_rows=len(phi),density_fractions=len(f),stable_samples=int((rates<-1e-7).sum()),
                 unstable_samples=int((rates>1e-7).sum()),marginal_samples=int((abs(rates)<=1e-7).sum()),
                 boundary_points=curves,precision_checks=selected)
    (DEST/'compensated_stability_map.json').write_text(json.dumps(summary,indent=2))
    plot(phi,phase_edges,f_edges,rates,curves);report(summary)

def plot(phi,phase_edges,f_edges,rates,curves):
    a,b=bounds(phase_edges);n=a[:,None]+(b-a)[:,None]*f_edges
    X=drive(n,phase_edges[:,None]);Y=np.broadcast_to(phase_edges[:,None],X.shape)
    labels=np.where(rates<-1e-7,0,np.where(rates>1e-7,1,2))
    colors=['#26866d','#cf6b4b','#edc64b'];cmap=ListedColormap(colors);norm=BoundaryNorm([-.5,.5,1.5,2.5],3)
    handles=[Patch(facecolor=colors[0],label='Stable compensated intermediate state'),Patch(facecolor=colors[1],label='Unstable intermediate state'),Patch(facecolor='#e7e7e7',label='No separate intermediate root'),Patch(facecolor=colors[2],label='Marginal / unresolved')]
    # Refined zero-growth points anchor interpolated boundaries and avoid the
    # sawtooth artifacts of coloring the nonuniform parameter mesh cellwise.
    centers=(f_edges[:-1]+f_edges[1:])/2
    source=np.load(DEST/'compensated_stability_map.npz')
    eps=source['epsilon'];xp=eps.ravel();yp=np.repeat(phi,rates.shape[1]);zp=rates.ravel()
    if curves:
        xp=np.r_[xp,[r['epsilon'] for r in curves]];yp=np.r_[yp,[r['phi'] for r in curves]];zp=np.r_[zp,np.zeros(len(curves))]
    triangles=Triangulation(xp,yp)
    contour_cmap=ListedColormap([colors[0],colors[2],colors[1]])
    densephi=np.linspace(0,CUSP,2001);lo,hi=bounds(densephi)
    lowdrive=drive(hi,densephi);highdrive=drive(lo,densephi)
    for name,limits in PLOT_RANGES:
        fig,ax=plt.subplots(figsize=(11,7),layout='constrained');ax.set_facecolor('#e7e7e7')
        ax.tricontourf(triangles,zp,levels=[-10,-1e-7,1e-7,10],cmap=contour_cmap)
        ax.fill_betweenx(densephi,0,lowdrive,color='#e7e7e7',zorder=3)
        ax.fill_betweenx(densephi,highdrive,85,color='#e7e7e7',zorder=3)
        ax.plot(X[:,0],phase_edges,color='k',lw=.7,zorder=4);ax.plot(X[:,-1],phase_edges,color='k',lw=.7,zorder=4)
        if curves:ax.scatter([r['epsilon'] for r in curves],[r['phi'] for r in curves],s=1,color='#25343d',alpha=.6,zorder=4)
        ax.axvline(38.51,color='#25343d',ls='--',lw=.9,zorder=4)
        ax.set(xlim=limits[:2],ylim=limits[2:],xlabel=r'Drive $\epsilon/J$',ylabel=r'Hopping phase $\phi$',title=PLOT_TITLE or f'Local stability of the compensated homogeneous intermediate state | OBC, N={N}')
        ticks=np.arange(0,7)*np.pi/6;labels_tick=['0',r'$\pi/6$',r'$\pi/3$',r'$\pi/2$',r'$2\pi/3$',r'$5\pi/6$',r'$\pi$']
        ax.set_yticks(ticks,labels_tick);ax.set_ylim(limits[2:]);ax.legend(handles=handles,loc='upper right',fontsize=9)
        ax.grid(alpha=.13);fig.savefig(f'figures/{FIGURE_PREFIX}_{name}.png',dpi=200)
        fig.savefig(f'figures/{FIGURE_PREFIX}_{name}.pdf')

def report(summary):
    text=f'''## Two-dimensional map: compensation-stabilizable intermediate solutions

The requested map uses **epsilon on the horizontal axis and phi on the vertical axis**, retaining N=50, J=1, U=-0.0002, V=0, Delta=0.5 and kappa=1. Green means the intermediate bulk root exists and its **homogeneous OBC BdG operator is strictly stable**. The exact fixed endpoint compensation therefore constructs a locally stable homogeneous mean-field BHM solution there. Orange means an existing intermediate target is unstable; grey means no separate intermediate root; yellow marks growth magnitudes <=1e-7 J, numerically marginal/unresolved.

![Full epsilon-phi map of compensation-stabilizable intermediate mean-field states, N=50 OBC](../figures/compensated_stability_map_full.png)

![Zoom over epsilon=30-40, with the epsilon=38.51 slice indicated by a dashed vertical line](../figures/compensated_stability_map_zoom.png)

This is **not a map of successful adiabatic preparation or robustness to imperfect compensation**. Green establishes exact stationarity and local BdG stability with root-specific fixed drives. Vacuum/ramp reachability can be more restrictive, as the separate simulations above show. It does not imply three stable roots for one fixed endpoint-drive setting.

### Existence region and sampling

The intermediate density lies between the positive turning densities of the bulk cubic. With c(phi)=0.4998+2 cos(phi), u=-0.0004 and gamma=0.5, these solve 3u^2 n^2+4cu n+c^2+gamma^2=0. Positive-density bistability requires c>sqrt(3)gamma, restricting the existence lobe to phi<{summary['cusp_degrees']:.6f} degrees. The lobe starts at the cusp drive epsilon=21.93456688 and reaches upper drive epsilon=79.51348142 at phi=0, so the full axes 0<=epsilon<=85 and 0<=phi<=pi cover all possible three-positive-root solutions for these parameters.

We parameterize the intermediate branch directly by its density between the two turning points. There are {summary['phase_rows']} phase rows, uniformly spaced over the existence sector (spacing about 0.2483 degrees), and {summary['density_fractions']} nonuniform density fractions per row, with extra samples near both root mergers. Each density is converted to epsilon through the stationary cubic and gets its own OBC spectrum. This curved sampling avoids missing thin stable strips next to an existence curve, which a coarse rectangular epsilon grid could miss. The shading interpolates growth rates on a triangulation augmented with refined zero-growth points; the analytical existence curves mask the outside region. Black outer curves mark root-existence limits, dark internal points mark refined zero-growth boundaries, and the dashed vertical line is epsilon=38.51.

The calculation has {summary['samples']} stationary-background spectra: {summary['stable_samples']} stable, {summary['unstable_samples']} unstable and {summary['marginal_samples']} marginal samples. These sample counts are not physical parameter-space areas because the density grid is nonuniform. Root-existence boundaries are excluded as distinct intermediate solutions because roots coincide there. The phi=0 boundary is approached by the mesh and checked separately at selected densities.

### Stability and verification

At each point d=Delta+U+4Un and |g|=|2U|n define the uniform OBC operator. A global Nambu rotation removes the constant pairing phase without changing the spectrum. Modes evolve as exp(-i omega t), so max Im(omega)<0 is stable. Endpoint compensation is calculated from the complex homogeneous alpha0, not optimized: pump_1=epsilon+J exp(+i phi) alpha0, pump_N=epsilon+J exp(-i phi) alpha0. Fixed additive compensation contributes no BdG derivatives. No Gaussian covariance equations or feedback are used.

All zero-growth sign changes along each sampled density row are refined with a scalar root solve. Representative points at phi=0,20,40,60,75 degrees, with densities near both turning boundaries and in the middle, are independently checked by characteristic-polynomial/Routh counts at 100 and 160 digits. Their stable/unstable labels agree with both precisions, and the corresponding complex compensated mean-field residuals are below 1e-8. The earlier high-precision epsilon=38.51 phase-slice checks remain independent validation. Selected checks do not certify every colored cell; near boundaries the plot is a finite-resolution numerical map, not a rigorous continuum certificate. N is fixed at 50 and no infinite-chain stability conclusion is claimed.

Reproduce with `python src/compensated_stability_map.py`. Curved-grid densities, drives, phases and growth rates are saved in data/generated/compensation_n50/compensated_stability_map.npz; counts, refined boundary points, parameters and precision checks are in compensated_stability_map.json. Full and zoom figures are provided as PNG and vector PDF files under figures/.
'''
    p=Path('docs/uniform_topological_segment_stability.md');base=p.read_text(encoding='utf-8').split('## Two-dimensional map: compensation-stabilizable intermediate solutions')[0].rstrip();p.write_text(base+'\n\n'+text,encoding='utf-8')

if __name__=='__main__':main()
