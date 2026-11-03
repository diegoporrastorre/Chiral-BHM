"""Reproduce fixed-background phase sweep and figures; no stationary-state solve."""
from pathlib import Path
import json
import numpy as np
from scipy.linalg import eig, eigvals
from scipy.optimize import brentq
from bulk_stability import matrix, bulk_threshold, local_response

ROOT=Path(__file__).resolve().parents[1]


def run():
    rows=[]
    for deg in sorted(set(list(range(0,91,5))+[50,52,54,56,58,60,62,64])):
        phi=np.deg2rad(deg)
        args=dict(J=1.,phi=phi,delta=.5,g=.4)
        threshold,candidates=bulk_threshold(**args,steps=1000)
        refined,_=bulk_threshold(**args,steps=2400,epsilon=1e-8)
        dmin=max(0,.5-2*abs(np.cos(phi)))
        pbc=float(np.sqrt(max(0,.4**2-dmin**2)))
        row=dict(degrees=deg,bulk=threshold,pbc=pbc,
                 pinch_refinement_difference=abs(threshold-refined),candidates=candidates,obc={})
        for N in [40,80,120,160]:
            H=matrix(N,**args)
            w,left,right=eig(H,left=True,right=True)
            idx=np.argmax(w.imag)
            cond=1/abs(np.vdot(left[:,idx],right[:,idx]))
            row['obc'][str(N)]=dict(threshold=float(w[idx].imag),condition=float(cond),
                residual=float(np.linalg.norm(H@right[:,idx]-w[idx]*right[:,idx])/np.linalg.norm(H)))
        rows.append(row)
        print(deg,'bulk',round(threshold,8),'obc160',round(row['obc']['160']['threshold'],8),flush=True)
    bulk_cross=brentq(lambda p:bulk_threshold(1,p,.5,.4,steps=1000)[0]-.35,np.deg2rad(30),np.deg2rad(70),xtol=1e-9)
    crossings={}
    for N in [40,80,120]:
        crossings[str(N)]=float(np.rad2deg(brentq(lambda p:eigvals(matrix(N,1,p,.5,.4)).imag.max()-.35,
                  np.deg2rad(30),np.deg2rad(70),xtol=1e-8)))
    impulse=[]
    for deg in [20,40,60,70,80]:
        ts=np.linspace(15,45,601)
        values=local_response(ts,1,np.deg2rad(deg),.5,.4,nk=8192)
        values2=local_response(ts,1,np.deg2rad(deg),.5,.4,nk=16384)
        fits=[]
        for lo,hi in [(15,35),(25,45),(15,45)]:
            mask=(ts>=lo)&(ts<=hi)
            fits.append(float(np.polyfit(ts[mask],np.log(values[mask]*np.sqrt(ts[mask])),1)[0]))
        impulse.append(dict(degrees=deg,fits=fits,quadrature_change=float(np.max(abs(values-values2)/values2))))
    result=dict(parameters=dict(J=1,delta=.5,g=.4,gamma=.35,kappa=.7),
                rows=rows,bulk_crossing_degrees=float(np.rad2deg(bulk_cross)),
                obc_crossing_degrees=crossings,impulse=impulse,
                stability_tolerance=1e-8,numpy_version=np.__version__)
    folder=ROOT/'data/generated';folder.mkdir(parents=True,exist_ok=True)
    (folder/'bdg_phi_sweep.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print('CROSSINGS',result['bulk_crossing_degrees'],crossings,flush=True)
    print('IMPULSE',impulse,flush=True)
    return result


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows=result['rows'];x=np.array([r['degrees'] for r in rows])
    fig,axes=plt.subplots(2,1,figsize=(8,7.5),layout='constrained',sharex=True)
    for ax,gamma in zip(axes,[.35,.1]):
        ax.plot(x,[r['pbc']-gamma for r in rows],color='0.35',ls='--',label='PBC, continuous k')
        ax.plot(x,[r['bulk']-gamma for r in rows],color='black',lw=2,label='Bulk saddle prediction')
        for N,marker in [(40,'o'),(80,'s'),(160,'^')]:
            vals=np.array([r['obc'][str(N)]['threshold']-gamma for r in rows])
            sensitive=np.array([r['obc'][str(N)]['condition']>1e9 for r in rows])
            line,=ax.plot(x,np.where(sensitive,np.nan,vals),marker=marker,ms=4,lw=.7,label=f'OBC, N={N}')
            ax.scatter(x[sensitive],vals[sensitive],marker='x',s=22,color=line.get_color(),alpha=.55)
        ax.axhline(0,color='0.5',lw=.8)
        ax.set(ylabel=r'Growth rate / $J$',title=fr'$\kappa/J={2*gamma:g}$',xlim=(0,90))
        ax.grid(alpha=.2)
    axes[0].legend(fontsize=8)
    axes[1].set_xlabel(r'Hopping phase $\phi$ (degrees)')
    fig.suptitle(r'Fixed background: $\tilde\Delta/J=0.5$, $|g_s|/J=0.4$, $g_c=0$')
    fig.supxlabel('Crosses: eigenvalue condition number > 1e9; interpret cautiously',fontsize=8)
    folder=ROOT/'figures';folder.mkdir(exist_ok=True)
    for suffix in ['png','pdf']:
        fig.savefig(folder/f'bdg_phi_stability.{suffix}',dpi=180)
    plt.close(fig)


def boundary_counterexample():
    """Independent propagation check where the bulk criterion fails for OBC."""
    from scipy.sparse import csr_matrix
    from scipy.sparse.linalg import expm_multiply
    args=dict(J=1.,phi=np.deg2rad(80),delta=.5,g=.4,gamma=.1)
    rows=[]
    for N in [40,80,120]:
        H=matrix(N,**args)
        w,left,right=eig(H,left=True,right=True)
        ix=np.argmax(w.imag)
        initial=np.zeros(2*N,dtype=complex)
        initial[N//2]=initial[N+N//2]=1/np.sqrt(2)
        ts=np.linspace(0,250,501)
        evolved=expm_multiply(csr_matrix(-1j*H),initial,start=0,stop=250,num=501)
        norms=np.linalg.norm(evolved,axis=1)
        fit=float(np.polyfit(ts[ts>=150],np.log(norms[ts>=150]),1)[0])
        # Explicit tiny deterministic perturbation checks spectral sensitivity.
        perturbed=H.copy();perturbed[0,0]+=1e-10
        sensitivity=abs(eigvals(perturbed).imag.max()-w[ix].imag)
        rows.append(dict(N=N,growth=float(w[ix].imag),omega_real=float(w[ix].real),
                         condition=float(1/abs(np.vdot(left[:,ix],right[:,ix]))),
                         time_fit_150_250=fit,perturbation_shift=float(sensitivity)))
    bulk,_=bulk_threshold(1,args['phi'],.5,.4)
    result=dict(parameters=args,bulk_growth=bulk-.1,rows=rows)
    (ROOT/'data/generated/bdg_boundary_counterexample.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2),flush=True)
    return result


if __name__=='__main__':
    plot(run())
    boundary_counterexample()
