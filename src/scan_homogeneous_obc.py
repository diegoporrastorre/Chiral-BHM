"""N=100 homogeneous mean-field scan; run from project root."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigvals, svd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mean_field import Chain

OUT = Path('data/generated/homogeneous_obc_n100')
OUT.mkdir(parents=True, exist_ok=True)

def evolve(m, a, duration=400, rtol=2e-8, atol=2e-9):
    def rhs(t, a):
        f = m.pump + (m.delta+m.U-.5j*m.kappa+2*m.U*abs(a)**2)*a
        f[:-1] += np.exp(-1j*m.phi)*a[1:]
        f[1:] += np.exp(1j*m.phi)*a[:-1]
        return -1j*f
    sol = solve_ivp(rhs, [0,duration], a, method='DOP853',rtol=rtol,atol=atol,
                    t_eval=np.linspace(duration-50,duration,101))
    return sol.y[:,-1], float(np.max(np.std(abs(sol.y)**2,axis=1))), sol

def diagnose(m, a):
    n=abs(a)**2; d=m.delta+m.U+4*m.U*n; g=2*m.U*a*a
    x0=abs(g)**2-.25-(d+1)**2; xp=abs(g)**2-.25-(d-1)**2
    nu=np.where(x0*xp<0,np.sign(d).astype(int),0)
    # Slow variation in BOTH complex pairing and detuning, measured in units J.
    variation=np.maximum(abs(np.diff(d)),abs(np.diff(g)))
    eta=np.maximum(np.r_[0,variation],np.r_[variation,0])
    L=m.bogoliubov(a)
    w=eigvals(L)
    u,s,vh=svd(-L)
    v=vh.conj().T[:,-1]; uu=u[:,-1]
    return dict(density=n.tolist(),alpha_real=a.real.tolist(),alpha_imag=a.imag.tolist(),
                residual=float(max(abs(m.force(a)))),growth=float(max(w.imag)),
                nu=nu.tolist(),variation=eta.tolist(),smin=float(s[-1]),snext=float(s[-2]),
                singular_centers=[float(np.dot(np.arange(1,m.N+1),abs(z[:m.N])**2+abs(z[m.N:])**2)) for z in (v,uu)])

def main():
    records=[]
    for direction, grid in [('up',np.arange(20,60.001,.25)),('down',np.arange(60,19.999,-.25))]:
        a=np.zeros(100,complex)
        for eps in grid:
            m=Chain(N=100,drive=float(eps))
            r=m.stationary(a,residual_tol=1e-8)
            stable=False
            if r['converged']:
                stable=max(eigvals(m.bogoliubov(r['alpha'])).imag)<-1e-7
            temporal_std=None
            if not stable:
                a,temporal_std,_=evolve(m,a)
                r=m.stationary(a,residual_tol=1e-8)
            if r['converged']:
                a=r['alpha']; rec=diagnose(m,a); rec['stationary']=True
            else:
                rec=dict(stationary=False,residual=r['residual'],density=(abs(a)**2).tolist())
            rec.update(epsilon=float(eps),direction=direction,temporal_std=temporal_std)
            records.append(rec)
            print(direction,eps,rec['stationary'],rec.get('growth'),flush=True)
        (OUT/'scan.json').write_text(json.dumps(records,indent=2))
    # Independent vacuum starts at representative drives; longer time and tighter tolerances.
    checks=[]
    for eps in [20,30,35,38,40,42,45,50,60]:
        m=Chain(N=100,drive=eps)
        a,std,sol=evolve(m,np.zeros(100,complex),duration=1200,rtol=2e-10,atol=2e-11)
        r=m.stationary(a,residual_tol=1e-8)
        rec=diagnose(m,r['alpha']) if r['converged'] else dict(residual=r['residual'])
        rec.update(epsilon=eps,stationary=r['converged'],temporal_std=std,
                   end_force=float(max(abs(m.force(a)))))
        checks.append(rec)
    (OUT/'checks.json').write_text(json.dumps(checks,indent=2))
    fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for col,direction in enumerate(['up','down']):
        rr=sorted([r for r in records if r['direction']==direction],key=lambda r:r['epsilon'])
        density=np.array([r['density'] if r['stationary'] else [np.nan]*100 for r in rr])
        nu=np.array([r.get('nu',[np.nan]*100) for r in rr])
        for row,z,label in [(0,density,'density'),(1,nu,'local winding at frequency zero')]:
            im=ax[row,col].imshow(z.T,origin='lower',aspect='auto',extent=[19.875,60.125,.5,100.5],cmap='viridis' if row==0 else 'coolwarm')
            fig.colorbar(im,ax=ax[row,col],label=label)
            ax[row,col].set(xlabel='drive epsilon / J',ylabel='site j',title=direction+' continuation')
    fig.savefig('figures/homogeneous_obc_n100_scan.png',dpi=180)
    fig,ax=plt.subplots(2,1,figsize=(10,7),layout='constrained')
    for eps in [20,30,35,38,40,42,45,50,60]:
        r=next(r for r in records if r['direction']=='up' and r['epsilon']==eps)
        if r['stationary']: ax[0].plot(range(1,101),r['density'],label=str(eps))
    ax[0].legend(ncol=5,title='epsilon / J');ax[0].set(xlabel='site j',ylabel='density')
    for direction in ['up','down']:
        rr=[r for r in records if r['direction']==direction and r['stationary']]
        ax[1].plot([r['epsilon'] for r in rr],[r['growth'] for r in rr],label=direction)
    ax[1].axhline(0,color='k',ls='--');ax[1].legend();ax[1].set(xlabel='epsilon / J',ylabel='max Im omega / J')
    fig.savefig('figures/homogeneous_obc_n100_profiles.png',dpi=180)

if __name__=='__main__': main()
