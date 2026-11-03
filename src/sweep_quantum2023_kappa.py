"""Exact one-way endpoint of Gomez-Leon et al., Quantum 7, 1016 (2023).

Uses their dynamical-matrix conventions, Eqs. (7)-(10), (28)-(31), P=0.
No stationary background is solved. kappa is the paper's gamma.
"""
from pathlib import Path
import json
import numpy as np
from scipy.linalg import eigvals, svdvals, expm, block_diag

ROOT = Path(__file__).resolve().parents[1]


def channels(n, kappa, pbc=False):
    # (S u)_j=u_{j+1}, so S has Bloch symbol exp(ik).
    s = np.diag(np.ones(n-1), 1)
    if pbc:
        s[-1, 0] = 1
    return 1j*((1-kappa/2)*np.eye(n)+2*s), -1j*((1+kappa/2)*np.eye(n)+2*s.T)


def nambu_matrix(n, kappa):
    s = np.diag(np.ones(n-1), 1)
    a = 1j*(s-s.T)  # phi=pi/2 represented exactly, not exp(i*pi/2)
    b = np.eye(n)+s+s.T
    return np.block([[a-1j*kappa/2*np.eye(n), b], [-b, a-1j*kappa/2*np.eye(n)]])


def channel_rotation(n):
    eye = np.eye(n)
    return np.block([[eye, eye], [1j*eye, -1j*eye]])/np.sqrt(2)


def winding(kappa, nk):
    if min(abs(kappa-2), abs(kappa-6)) < 1e-12:
        return None
    k = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    plus = 1j*(1-kappa/2)+2j*np.exp(1j*k)
    minus = -1j*(1+kappa/2)-2j*np.exp(-1j*k)
    q = plus*minus
    return int(np.rint(np.angle(np.roll(q,-1)/q).sum()/(2*np.pi)))


def propagator(n, kappa, t):
    # Exact finite Taylor polynomial of each nilpotent shift.
    s = np.diag(np.ones(n-1), 1)
    term = np.eye(n)
    pp = term.copy()
    pm = term.copy()
    for d in range(1, n):
        term = term @ s * (2*t/d)
        pp += term
        pm += (-1)**d * term.T
    return block_diag(np.exp((1-kappa/2)*t)*pp, np.exp((-1-kappa/2)*t)*pm)


def run():
    rows = []
    for kap in np.linspace(0,8,81):
        w = winding(kap,4096)
        assert w == winding(kap,8192)
        assert w == (None if kap in (2,6) else int(2<kap<6))
        rows.append(dict(kappa=float(kap), bulk_growth=1-kap/2,
                         obc_growth=1-kap/2, pbc_growth=3-kap/2, winding=w))
    checks=[]
    for n in [8,16,32,64,96]:
        kap=3.
        hp,hm=channels(n,kap)
        hc=block_diag(hp,hm)
        h=nambu_matrix(n,kap)
        u=channel_rotation(n)
        residual=float(np.linalg.norm(u.conj().T@h@u-hc,ord=np.inf))
        assert residual<1e-13
        rate=float(eigvals(hc).imag.max())
        assert abs(rate+.5)<1e-12
        checks.append(dict(N=n,rotation_residual=residual,channel_eig_growth=rate,
                           raw_nambu_eig_growth=float(eigvals(h).imag.max())))
    svds=[]
    for kap in [1.,2.,3.,4.,5.,5.5,6.,7.]:
        for n in [8,16,32,64,96]:
            hp,hm=channels(n,kap)
            vals=np.sort(np.r_[svdvals(hp),svdvals(hm)])
            svds.append(dict(kappa=kap,N=n,smin=float(vals[0]),snext=float(vals[1])))
    evolution=[]
    for kap in [1.,2.,3.,5.,7.]:
        n=8
        h=nambu_matrix(n,kap)
        u=channel_rotation(n)
        for t in [0.,1.,5.,10.,20.,40.,80.]:
            exact=propagator(n,kap,t)
            numerical=u.conj().T@expm(-1j*h*t)@u
            err=float(np.linalg.norm(numerical-exact)/max(1.,np.linalg.norm(exact)))
            # Long-time full-basis expm also suffers from nonnormal roundoff.
            if t<=10:
                assert err<1e-7,(kap,t,err)
            evolution.append(dict(kappa=kap,N=n,t=t,norm=float(np.linalg.norm(exact,2)),
                                  expm_relative_error=err))
    pbc=[]
    for n in [16,32,64]:
        for kap in [1.,3.,5.,7.]:
            rate=float(max(eigvals(c).imag.max() for c in channels(n,kap,True)))
            assert abs(rate-(3-kap/2))<1e-12
            pbc.append(dict(N=n,kappa=kap,numerical_growth=rate))
    result=dict(parameters=dict(gs=1,gc=1,Delta=0,J=1,phi='pi/2',P=0),
                loss='kappa=paper gamma; amplitude damping kappa/2', reference_frequency=0,
                rows=rows,checks=checks,svds=svds,evolution=evolution,pbc=pbc)
    path=ROOT/'data/generated/quantum2023_kappa.json'
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2))
    print(json.dumps(dict(checks=checks,max_short_time_expm_error=max(r['expm_relative_error'] for r in evolution if r['t']<=10)),indent=2))
    return result


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(2,2,figsize=(11,7),layout='constrained')
    ax=axs[0,0]
    x=np.array([r['kappa'] for r in result['rows']])
    ax.axvspan(2,6,color='gold',alpha=.2,label='Nontrivial winding at frequency 0')
    ax.plot(x,1-x/2,label='Exact bulk = OBC (all finite N)')
    ax.plot(x,3-x/2,label='PBC')
    ax.axhline(0,color='k',lw=.7)
    ax.set(xlabel=r'$\kappa/J$',ylabel='Asymptotic growth rate / J')
    ax.legend(fontsize=8)
    ax=axs[0,1]
    for kap in [1.,3.,4.,5.,5.5,6.,7.]:
        rr=[r for r in result['svds'] if r['kappa']==kap]
        ax.semilogy([r['N'] for r in rr],[r['smin'] if r['smin']>0 else np.nan for r in rr],'o-',label=str(kap))
    ax.set(xlabel='N',ylabel='Smallest singular value at frequency 0')
    ax.legend(title=r'$\kappa/J$',fontsize=8,ncol=2)
    ax=axs[1,0]
    for kap in [1.,2.,3.,5.,7.]:
        rr=[r for r in result['evolution'] if r['kappa']==kap]
        ax.semilogy([r['t'] for r in rr],[r['norm'] for r in rr],'o-',label=str(kap))
    ax.set(xlabel='Jt',ylabel=r'$\|G^{\rm OBC}(t)\|_2$, N=8')
    ax.legend(title=r'$\kappa/J$',fontsize=8)
    ax=axs[1,1]
    rr=result['checks']
    ax.plot([r['N'] for r in rr],[r['raw_nambu_eig_growth'] for r in rr],'o-',label='Unstructured double precision eig')
    ax.axhline(-.5,label='Exact / channel-preserving eig',color='k')
    ax.axhline(0,lw=.7,color='gray')
    ax.set(xlabel='N',ylabel=r'Computed growth at $\kappa=3J$')
    ax.legend(fontsize=8)
    fig.suptitle(r'Quantum 2023 BdG convention: $g_s=g_c=J=1$, $\Delta=0$, $\phi=\pi/2$, $P=0$')
    fig.savefig(ROOT/'figures/quantum2023_kappa.png',dpi=180)
    plt.close(fig)


if __name__=='__main__':
    plot(run())
