"""Rassaert-inspired boundary envelopes for the fixed Quantum-2023 BdG chain.

This is a prescribed quadratic taper, not the paper's nonlinear solution.
Bulk: gs=gc=J=1, Delta=.1, phi=pi/2, amplitude decay kappa/2.
"""
from pathlib import Path
import json
import numpy as np
from scipy.linalg import eig, eigvals, svdvals
from scipy.sparse.linalg import expm_multiply
from sweep_quantum2023_kappa import channel_rotation
from sweep_quantum2023_detuned import matrix as sharp_matrix, recurrence_root

ROOT=Path(__file__).resolve().parents[1]


def profiles(n,width,kind='paired'):
    d=np.minimum(np.arange(1,n+1),np.arange(n,0,-1))
    s=np.ones(n) if width==0 else np.tanh(d/width)
    delta=.1*s*s
    gs=s*s if kind=='paired' else np.ones(n)
    gc=s[:-1]*s[1:] if kind=='paired' else np.ones(n-1)
    return s,delta,gs,gc


def matrix(n,width,kappa=0,kind='paired'):
    s,delta,gs,gc=profiles(n,width,kind)
    shift=np.diag(np.ones(n-1),1)
    hopping=1j*(shift-shift.T)
    a=hopping+np.diag(delta)
    b=np.diag(gs)+np.diag(gc,1)+np.diag(gc,-1)
    return np.block([[a,b],[-b,-a.conj()]])-1j*kappa/2*np.eye(2*n)


def spectrum(n,width,kind):
    h=matrix(n,width,kind=kind)
    w,vl,vr=eig(h,left=True,right=True)
    j=np.argmax(w.imag)
    u=channel_rotation(n)
    wc=eigvals(u.conj().T@h@u)
    norms=np.linalg.norm(vl,axis=0)*np.linalg.norm(vr,axis=0)
    cond=norms/np.abs(np.sum(vl.conj()*vr,axis=0))
    # First-order sensitivity indicator, not a rigorous error bound.
    sensitivity=float(np.finfo(float).eps*np.linalg.norm(h,2)*cond[j])
    s,delta,gs,gc=profiles(n,width,kind)
    return dict(N=n,width=width,kind=kind,growth=float(w[j].imag),
                critical_kappa=float(2*w[j].imag),frequency=float(w[j].real),
                dominant_condition=float(cond[j]),sensitivity_indicator=sensitivity,
                basis_difference=float(abs(w[j].imag-wc.imag.max())),
                eigen_residual=float(np.linalg.norm(h@vr[:,j]-w[j]*vr[:,j])/np.linalg.norm(h)),
                center_gs=float(gs[n//2]),center_delta=float(delta[n//2]),
                bulk_like_sites=int(np.sum(s*s>.99)))


def run():
    rows=[]
    for n in [16,24,32,40,64,80]:
        assert np.max(abs(matrix(n,0)-sharp_matrix(n)))<1e-14
        for kind in ['paired','detuning_only']:
            for width in [0,1,2,4,6,8]:
                r=spectrum(n,width,kind)
                if width==0:
                    r['recurrence_growth']=float(recurrence_root(n))
                rows.append(r)
                print(n,kind,width,round(r['critical_kappa'],7),r['sensitivity_indicator'],flush=True)
    result=dict(parameters=dict(gs=1,gc=1,J=1,Delta=.1,phi='pi/2',P=0),
                loss='uniform kappa; amplitude damping kappa/2',rows=rows)
    (ROOT/'data/generated/bdg_damped_boundaries.json').write_text(json.dumps(result,indent=2))
    return result


def refine(result):
    from banded_stability_decimal import threshold
    checks=[]
    cases=[(n,w,kind) for n in [16,32,40,64,80] for w in [2,4,6,8]
           for kind in ['paired','detuning_only']]
    cases += [(n,w,'paired') for n in [40,64,80] for w in [12,16]]
    for n,w,kind in cases:
        coarse=threshold(n,w,kind,120)
        fine=threshold(n,w,kind,200)
        assert abs(coarse['critical_kappa']-fine['critical_kappa'])<2e-8
        s,delta,gs,gc=profiles(n,w,kind)
        fine.update(precision_comparison=abs(coarse['critical_kappa']-fine['critical_kappa']),
                    center_gs=float(gs[n//2]),center_delta=float(delta[n//2]),
                    bulk_like_sites=int(np.sum(s*s>.99)))
        checks.append(fine)
        print('REFINED',n,w,kind,fine['critical_kappa'],flush=True)
        result['refined']=checks
        (ROOT/'data/generated/bdg_damped_boundaries.json').write_text(json.dumps(result,indent=2))
    return result


def dynamics(result,save=True):
    checks=[]; n=40; kap=3.5
    for kind,width in [('paired',0),('paired',4),('paired',6),('paired',8),('detuning_only',8)]:
        h=matrix(n,width,kap,kind)
        init=np.zeros(2*n,dtype=complex)
        init[n//2]=init[n+n//2]=1/np.sqrt(2)
        expected=(float(recurrence_root(n)) if width==0 else
                  next(r['critical_kappa']/2 for r in result['refined']
                       if r['N']==n and r['width']==width and r['kind']==kind))-kap/2
        # Remove a known scalar growth factor to prevent overflow/underflow.
        t=np.linspace(0,2000,1001)
        y=expm_multiply(-1j*h-expected*np.eye(2*n),init,start=0,stop=2000,num=1001)
        lognorms=np.log(np.linalg.norm(y,axis=1))+expected*t
        rate=float(np.polyfit(t[-501:],lognorms[-501:],1)[0])
        vals=svdvals(h)
        checks.append(dict(N=n,kappa=kap,kind=kind,width=width,expected=expected,fit=rate,
                           times=t.tolist(),log_norms=lognorms.tolist(),smin=float(vals[-1]),snext=float(vals[-2])))
        print('DYNAMICS',kind,width,rate,expected,flush=True)
        assert abs(rate-expected)<2e-3
    result['dynamics']=checks
    if save:
        (ROOT/'data/generated/bdg_damped_boundaries.json').write_text(json.dumps(result,indent=2))
    return result


def validate_characteristic():
    from decimal import Decimal,localcontext
    from banded_stability_decimal import characteristic_banded,envelope_generator
    records=[]
    with localcontext() as ctx:
        ctx.prec=80
        for kind in ['paired','detuning_only']:
            n=8; width=3
            a=envelope_generator(n,width,kind)
            af=np.array([[float(v) for v in row] for row in a])
            wa=eigvals(af); wh=-1j*eigvals(matrix(n,width,kind=kind))
            err=float(max(min(abs(wh-z)) for z in wa))
            assert err<1e-9
            coeff=characteristic_banded(a)
            x=Decimal('3.7'); value=Decimal(0)
            for c in coeff[::-1]: value=value*x+c
            dense=float(np.linalg.det(3.7*np.eye(2*n)-af))
            det_error=abs(dense-float(value))/abs(float(value))
            assert det_error<1e-10
            records.append(dict(kind=kind,spectrum_error=err,det_relative_error=det_error))
    return records


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(2,2,figsize=(11,7),layout='constrained')
    ax=axs[0,0]
    for width in [2,4,6,8]:
        s,d,gs,gc=profiles(40,width)
        ax.plot(np.arange(1,41),gs,label=f'N0={width}')
    ax.set(xlabel='Site j',ylabel=r'$g_{s,j}/J=\Delta_j/(0.1J)$',title='Prescribed boundary envelope, N=40')
    ax.legend(fontsize=8)
    ax=axs[0,1]
    for n in [16,32,40,64,80]:
        rr=sorted([r for r in result['refined'] if r['N']==n and r['kind']=='paired'],key=lambda r:r['width'])
        xs=[0]+[r['width'] for r in rr]
        ys=[2*float(recurrence_root(n))]+[r['critical_kappa'] for r in rr]
        ax.plot(xs,ys,'o-',label=f'N={n}')
    ax.axhline(2,color='gray',ls='--',label='Original homogeneous bulk threshold')
    ax.set(xlabel='Boundary width N0',ylabel=r'Critical $\kappa/J$',title='Pairing + detuning taper (high precision)')
    ax.legend(fontsize=7,ncol=2)
    ax=axs[1,0]; x=np.linspace(1.5,5,100)
    ax.plot(x,float(recurrence_root(40))-x/2,'k-',label='Sharp OBC')
    for kind,width in [('paired',4),('paired',6),('paired',8),('detuning_only',8)]:
        row=next(r for r in result['refined'] if r['N']==40 and r['width']==width and r['kind']==kind)
        ax.plot(x,row['critical_kappa']/2-x/2,label=f'{kind.replace("_"," ")}, N0={width}')
    ax.plot(x,1-x/2,'k--',label='Original local bulk')
    ax.axhline(0,color='gray',lw=.7)
    ax.set(xlabel=r'$\kappa/J$',ylabel='Growth rate / J',title='N=40: stability versus loss')
    ax.legend(fontsize=7)
    ax=axs[1,1]
    for row in result['dynamics']:
        times=np.array(row['times']); select=times<=600
        ax.semilogy(times[select],np.exp(np.array(row['log_norms'])[select]),label=f"{row['kind'].replace('_',' ')}, N0={row['width']}")
    ax.set(xlabel='Jt',ylabel='Physical impulse norm',title=r'$N=40$, $\kappa=3.5J$: direct dynamics')
    ax.legend(fontsize=7)
    fig.suptitle(r'Rassaert-inspired BdG taper: bulk $\Delta=0.1J$, $g_s=g_c=J=1$, $\phi=\pi/2$')
    fig.savefig(ROOT/'figures/bdg_damped_boundaries.png',dpi=180)
    plt.close(fig)


if __name__=='__main__':
    result=dynamics(refine(run()))
    result['characteristic_validation']=validate_characteristic()
    (ROOT/'data/generated/bdg_damped_boundaries.json').write_text(json.dumps(result,indent=2))
    plot(result)
