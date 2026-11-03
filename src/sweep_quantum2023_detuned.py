"""Detuned Quantum-2023 endpoint: fixed gs=gc=J=1, phi=pi/2, P=0."""
from pathlib import Path
from decimal import Decimal, localcontext
import json
import numpy as np
from scipy.linalg import eig, eigvals, svdvals
from scipy.sparse.linalg import expm_multiply
from sweep_quantum2023_kappa import nambu_matrix, channel_rotation, channels

ROOT=Path(__file__).resolve().parents[1]
DELTA=.1


def matrix(n,kappa=0,delta=DELTA):
    return nambu_matrix(n,kappa)+np.diag(np.r_[np.full(n,delta),np.full(n,-delta)])


def recurrence_root(n, precision=70):
    """Positive real characteristic root for even N; Decimal bisection.

    Determinant recurrence is independently derived, not an eigensolver fit.
    Other roots' spectral ordering is checked numerically at moderate N.
    """
    assert n%2==0
    with localcontext() as ctx:
        ctx.prec=precision
        d=Decimal('0.1')
        def det(lam):
            x=lam*lam
            prev=Decimal(1); cur=x-1+d*d
            for _ in range(2,n+1):
                prev,cur=cur,(x-5+d*d)*cur+4*(x-1)*prev
            return cur
        lo=Decimal(1); hi=(5-d*d).sqrt()
        assert det(lo)<0 and det(hi)>0
        for _ in range(precision*3):
            mid=(lo+hi)/2
            if det(mid)>0: hi=mid
            else: lo=mid
        return str((lo+hi)/2)


def local_response(t,nk=8192):
    k=np.arange(nk)*2*np.pi/nk
    b=1+2*np.cos(k)
    r=np.sqrt((b*b-DELTA**2).astype(complex))
    f=np.exp(2j*np.sin(k)*t)
    si=np.sinh(r*t)/r
    co=np.cosh(r*t)
    block=np.array([[np.mean(f*(co+b*si)),np.mean(-1j*DELTA*f*si)],
                    [np.mean(-1j*DELTA*f*si),np.mean(f*(co-b*si))]])
    return block


def right_half_plane_count(n, shift, precision=160):
    """Routh count for det[(s+shift)I-A], using Decimal arithmetic.

    No near-zero first-column entries are patched: those require more
    precision or a shift off the spectrum. Positive row scaling is allowed.
    """
    with localcontext() as ctx:
        ctx.prec=precision
        zero=Decimal(0)
        p=[Decimal(1)]; q=[Decimal('-.99'),Decimal(1)]
        for _ in range(2,n+1):
            r=[zero]*(len(q)+1)
            for j,v in enumerate(q): r[j]-=Decimal('4.99')*v; r[j+1]+=v
            for j,v in enumerate(p): r[j]-=4*v; r[j+1]+=4*v
            p,q=q,r
        # Ascending coefficients in lambda, then Horner shift lambda=s+shift.
        coeff=[zero]*(2*n+1)
        for j,v in enumerate(q): coeff[2*j]=v
        shifted=[coeff[-1]]
        c=Decimal(str(shift))
        for v in coeff[-2::-1]:
            new=[zero]*(len(shifted)+1)
            for j,z in enumerate(shifted): new[j]+=c*z; new[j+1]+=z
            new[0]+=v; shifted=new
        desc=shifted[::-1]
        a=desc[::2]; b=desc[1::2]+[zero]
        signs=[a[0]>0,b[0]>0]
        for _ in range(2,2*n+1):
            assert b[0]!=0
            r=[(b[0]*a[j+1]-a[0]*b[j+1])/b[0] for j in range(len(a)-1)]+[zero]
            scale=max(abs(v) for v in r)
            assert scale!=0
            r=[v/scale for v in r]
            signs.append(r[0]>0)
            a,b=b,r
        return sum(a!=b for a,b in zip(signs,signs[1:]))


def winding(kappa,nk):
    k=np.arange(nk)*2*np.pi/nk
    # det H = f0^2 - Delta^2 + G^2
    q=(-2*np.sin(k)-1j*kappa/2)**2-DELTA**2+(1+2*np.cos(k))**2
    return int(np.rint(np.angle(np.roll(q,-1)/q).sum()/(2*np.pi)))


def run():
    sizes=[4,8,16,24,32,48,64,96,160]
    checks=[]
    for n in sizes:
        h=matrix(n)
        w,vl,vr=eig(h,left=True,right=True)
        ix=np.argmax(w.imag)
        hp,hm=channels(n,0)
        hc=np.block([[hp,DELTA*np.eye(n)],[DELTA*np.eye(n),hm]])
        wc=eigvals(hc)
        u=channel_rotation(n)
        ref=float(recurrence_root(n,70))
        assert abs(ref-float(recurrence_root(n,100)))<1e-14
        counts=[]
        for precision in [160,240]:
            pair=[right_half_plane_count(n,ref+eps,precision) for eps in [-1e-8,1e-8]]
            assert pair==[1,0],(n,precision,pair)
            counts.append(pair)
        err=float(abs(w[ix].imag-ref))
        if n<=24: assert err<1e-7
        checks.append(dict(N=n,recurrence_growth=ref,critical_kappa=2*ref,
            raw_nambu_growth=float(w[ix].imag),raw_channel_growth=float(wc.imag.max()),
            routh_counts_below_above=counts,condition=float(1/abs(np.vdot(vl[:,ix],vr[:,ix]))),
            eigen_residual=float(np.linalg.norm(h@vr[:,ix]-w[ix]*vr[:,ix])/np.linalg.norm(h)),
            rotation_error=float(np.linalg.norm(u.conj().T@h@u-hc,np.inf))))
    large=[dict(N=n,growth=float(recurrence_root(n))) for n in [256,512,1024,2048]]
    lower=2*np.sqrt(1-DELTA**2); upper=2*np.sqrt(9-DELTA**2)
    rows=[]
    for kap in np.linspace(0,8,161):
        wind=winding(kap,8192)
        assert wind==winding(kap,16384)==int(lower<kap<upper)
        rows.append(dict(kappa=float(kap),winding=wind,bulk_growth=1-kap/2,
                         pbc_growth=np.sqrt(9-DELTA**2)-kap/2))
    response=[]
    for t in [1.,2.,4.,6.,8.,10.,12.]:
        r=local_response(t); fine=local_response(t,16384)
        response.append(dict(t=t,plus_scaled_real=float((fine[0,0]*np.exp(-t)).real),
                             quadrature_scaled_difference=float(np.linalg.norm(r-fine)*np.exp(-t))))
    svds=[]
    for kap in [1.,2.5,3.,4.,4.6,5.,5.8,6.5]:
        for n in [8,16,24,32,48,64]:
            s=svdvals(matrix(n,kap))
            svds.append(dict(kappa=kap,N=n,smin=float(s[-1]),snext=float(s[-2])))
    time_checks=[]
    n=16
    for kap in [3.,4.,5.]:
        gen=-1j*matrix(n,kap)
        init=np.zeros(2*n,dtype=complex); init[n//2]=1
        times=np.linspace(0,500,501)
        y=expm_multiply(gen,init,start=0,stop=500,num=501)
        norms=np.linalg.norm(y,axis=1)
        rate=float(np.polyfit(times[-81:],np.log(norms[-81:]),1)[0])
        expected=float(recurrence_root(n))-kap/2
        assert abs(rate-expected)<1e-5,(kap,rate,expected)
        time_checks.append(dict(N=n,kappa=kap,fit_growth=rate,expected=expected,
                               times=times.tolist(),norms=norms.tolist()))
    out=dict(parameters=dict(gs=1,gc=1,J=1,Delta=.1,phi='pi/2',P=0),
             conventions='kappa is Lindblad loss; i dt Psi=H Psi; omega_ref=0',
             topology_interval=[lower,upper],bulk_threshold=2.,pbc_threshold=upper,
             obc_limit_threshold=2*np.sqrt(5-DELTA**2),checks=checks,large=large,
             rows=rows,response=response,svds=svds,time_checks=time_checks)
    (ROOT/'data/generated/quantum2023_detuned.json').write_text(json.dumps(out,indent=2))
    print(json.dumps({k:out[k] for k in ['topology_interval','obc_limit_threshold','checks','large','response']},indent=2))
    return out


def plot(out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(2,2,figsize=(11,7),layout='constrained')
    ax=axs[0,0]; x=np.linspace(0,8,161)
    ax.axvspan(*out['topology_interval'],color='gold',alpha=.2,label='Nontrivial winding at frequency 0')
    ax.plot(x,1-x/2,'k--',label='Local bulk')
    for n in [8,16,64,160]:
        r=next(r for r in out['checks'] if r['N']==n)
        ax.plot(x,r['recurrence_growth']-x/2,label=f'OBC N={n} (recurrence)')
    ax.plot(x,out['pbc_threshold']/2-x/2,label='PBC')
    ax.axhline(0,color='gray',lw=.7); ax.set(xlabel=r'$\kappa/J$',ylabel='Growth rate / J'); ax.legend(fontsize=7)
    ax=axs[0,1]
    for key,label in [('recurrence_growth','70-digit characteristic recurrence'),('raw_nambu_growth','Raw Nambu eig'),('raw_channel_growth','Raw channel eig')]:
        ax.plot([r['N'] for r in out['checks']],[r[key] for r in out['checks']],'o-',label=label)
    ax.axhline(out['obc_limit_threshold']/2,color='gray',ls='--',label='Large-N band edge')
    ax.set(xlabel='N',ylabel='Lossless OBC growth / J'); ax.legend(fontsize=7)
    ax=axs[1,0]
    for kap in [1.,2.5,4.,4.6,5.,5.8,6.5]:
        rr=[r for r in out['svds'] if r['kappa']==kap]
        ax.semilogy([r['N'] for r in rr],[r['smin'] for r in rr],'o-',label=str(kap))
    ax.axhline(1e-14,color='gray',ls=':',label='Double precision caution')
    ax.set(xlabel='N',ylabel='Smallest singular value at frequency 0'); ax.legend(title=r'$\kappa/J$',fontsize=7,ncol=2)
    ax=axs[1,1]
    for r in out['time_checks']:
        ax.semilogy(r['times'],r['norms'],label=f"kappa={r['kappa']}, fitted rate {r['fit_growth']:.4f}")
    ax.set(xlabel='Jt',ylabel='Impulse norm, N=16'); ax.legend(fontsize=8)
    fig.suptitle(r'$g_s=g_c=J=1$, $\Delta=0.1J$, $\phi=\pi/2$, $P=0$: coupling the channels')
    fig.savefig(ROOT/'figures/quantum2023_detuned.png',dpi=180)
    plt.close(fig)


if __name__=='__main__':
    plot(run())
