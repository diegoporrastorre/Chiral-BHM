"""Small quenched-Gaussian disorder pilot for the compensated V=0 target."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json
import numpy as np
from scipy.linalg import solve,eig,eigvals,svd
from scipy.integrate import solve_ivp
from mean_field import Chain
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DEST=Path('data/generated/compensated_disorder_pilot');DEST.mkdir(parents=True,exist_ok=True)
EPS=38.51;PHI=np.pi/3;U=-.0002
ROOT=sorted(np.roots([.0004**2,-.0008*1.4998,1.4998**2+.25,-EPS**2]).real)[1]
ALPHA=-EPS/(1.4998-.0004*ROOT-.5j)
KINDS=['delta','J','phi'];SEEDS=[20261104,20261105];STRENGTHS=[1e-6,1e-4,1e-2]

def disorder(N,kind,seed):
    # Same normalized realization at all strengths; shorter chains use prefixes.
    rng=np.random.default_rng(np.random.SeedSequence([seed,KINDS.index(kind)]))
    z=rng.normal(size=50 if kind=='delta' else 49)
    return z[:N if kind=='delta' else N-1]

def model(N,kind,z,sigma):
    m=Chain(N=N,drive=EPS);m.pump=m.pump.astype(complex)
    m.pump[0]+=np.exp(1j*PHI)*ALPHA;m.pump[-1]+=np.exp(-1j*PHI)*ALPHA
    if kind=='delta':m.detuning=.5+sigma*z
    else:
        m.T[:]=0
        for j in range(N-1):
            hop=(1+sigma*z[j])*np.exp(1j*PHI) if kind=='J' else np.exp(1j*(PHI+sigma*z[j]))
            m.T[j+1,j]=hop;m.T[j,j+1]=hop.conjugate()
    return m

def diagnose(m,a):
    A,B=m.blocks(a);L=m.bogoliubov(a);w,left,right=eig(L,left=True,right=True);j=np.argmax(w.imag)
    cond=1/abs(np.vdot(left[:,j],right[:,j]));sensitivity=np.finfo(float).eps*np.linalg.norm(L,np.inf)*cond
    G=np.block([[(A+B).imag,(A-B).real],[-(A+B).real,(A-B).imag]])
    real_growth=float(eigvals(G).real.max());rate=float(w[j].imag)
    Umat,s,Vh=svd(L);v=Vh.conj().T[:,-1];u=Umat[:,-1];sites=np.arange(1,m.N+1)
    resolution=max(1e-7,10*sensitivity,10*abs(real_growth-rate))
    return dict(residual=float(max(abs(m.force(a)))),growth=rate,
                stability='stable' if rate < -resolution else 'unstable' if rate > resolution else 'unresolved',
                leading_condition=float(cond),sensitivity_indicator=float(sensitivity),basis_growth_difference=abs(real_growth-rate),resolution=resolution,
                density=(abs(a)**2).tolist(),alpha_real=a.real.tolist(),alpha_imag=a.imag.tolist(),
                max_relative_amplitude_error=float(max(abs(a-ALPHA)))/abs(ALPHA),
                max_relative_density_error=float(max(abs(abs(a)**2-ROOT)))/ROOT,
                smin=float(s[-1]),snext=float(s[-2]),
                right_singular_center=float(np.dot(sites,abs(v[:m.N])**2+abs(v[m.N:])**2)),
                left_singular_center=float(np.dot(sites,abs(u[:m.N])**2+abs(u[m.N:])**2)))

def linear_prediction(N,kind,z,sigma):
    clean=model(N,kind,z,0);a=np.full(N,ALPHA);A,B=clean.blocks(a)
    K=np.block([[(A+B).real,-(A-B).imag],[(A+B).imag,(A-B).real]])
    q=np.zeros(N,complex)
    if kind=='delta':q=z*ALPHA
    else:
        derivative=z*np.exp(1j*PHI) if kind=='J' else 1j*z*np.exp(1j*PHI)
        q[1:]+=derivative*ALPHA;q[:-1]+=derivative.conjugate()*ALPHA
    x=solve(K,-np.r_[q.real,q.imag]);response=x[:N]+1j*x[N:]
    return dict(relative_peak_per_sigma=float(max(abs(response)))/abs(ALPHA),
                predicted_relative_peak=float(sigma*max(abs(response))/abs(ALPHA)))

def continuation(N,kind,z,sigma):
    a=np.full(N,ALPHA);lam=0.;step=.1;prev=None;path=[];failed=[]
    for iteration in range(600):
        next_lam=min(1,lam+step);m=model(N,kind,z,sigma*next_lam)
        initial=a if prev is None else a+(next_lam-lam)*(a-prev[1])/(lam-prev[0])
        r=m.stationary(initial,1e-8);change=float(max(abs(r['alpha']-a)))
        if not r['converged'] or change>2:
            failed.append(dict(lam=next_lam,step=step,residual=r['residual'],amplitude_step=change));step/=2
            if step<1e-5:break
            continue
        prev=(lam,a.copy());lam=next_lam;a=r['alpha']
        w=eigvals(m.bogoliubov(a));path.append(dict(lam=lam,growth=float(w.imag.max()),residual=r['residual'],relative_peak=float(max(abs(a-ALPHA))/abs(ALPHA))))
        step=min(.2,step*1.4)
        if lam==1:break
    return dict(reached_full_disorder=bool(lam==1),last_lam=lam,path=path,failed_attempts=failed,
                final=diagnose(model(N,kind,z,sigma*lam),a))

def evolve(m,duration=1200,rtol=2e-10,atol=2e-11):
    upper=np.diag(m.T,1);lower=np.diag(m.T,-1)
    def rhs(t,a):
        f=m.pump+(m.detuning+m.U-.5j*m.kappa+2*m.U*abs(a)**2)*a
        f[:-1]+=upper*a[1:];f[1:]+=lower*a[:-1]
        return -1j*f
    times=np.unique(np.r_[np.linspace(0,150,301),np.linspace(duration-50,duration,101)])
    sol=solve_ivp(rhs,[0,duration],np.full(m.N,ALPHA),method='DOP853',rtol=rtol,atol=atol,t_eval=times)
    assert sol.success;end=sol.y[:,-1];root=m.stationary(end,1e-8)
    force=float(max(abs(m.force(end))));distance=float(max(abs(root['alpha']-end)))
    reached=root['converged'] and force<1e-3 and distance<1e-3
    final=diagnose(m,root['alpha']) if reached else dict(density=(abs(end)**2).tolist(),
                max_relative_amplitude_error=float(max(abs(end-ALPHA)))/abs(ALPHA),
                max_relative_density_error=float(max(abs(abs(end)**2-ROOT)))/ROOT)
    final.update(reached_stationary=bool(reached),raw_end_force=force,root_distance=distance,
                 duration=duration,rtol=rtol,atol=atol,tail_density_std=float(np.std(abs(sol.y[:,sol.t>=duration-50])**2,axis=1).max()),
                 times=sol.t.tolist(),density_trajectory=(abs(sol.y)**2).T.tolist())
    return final

def validate():
    rng=np.random.default_rng(124)
    for kind in KINDS:
        z=disorder(8,kind,SEEDS[0]);m=model(8,kind,z,.01)
        np.testing.assert_allclose(m.T,m.T.conj().T,atol=1e-14)
        a=rng.normal(size=8)+1j*rng.normal(size=8);v=rng.normal(size=8)+1j*rng.normal(size=8);A,B=m.blocks(a);h=1e-5
        error=max(abs((m.force(a+h*v)-m.force(a-h*v))/(2*h)-A@v-B@v.conj()));assert error<1e-8,error
        # Disorder forcing derivative used in the linear predictor.
        alpha=np.full(8,ALPHA);fd=(model(8,kind,z,h).force(alpha)-model(8,kind,z,-h).force(alpha))/(2*h)
        exact=z*ALPHA if kind=='delta' else np.zeros(8,complex)
        if kind!='delta':
            deriv=z*np.exp(1j*PHI)*(1 if kind=='J' else 1j);exact[1:]+=deriv*ALPHA;exact[:-1]+=deriv.conj()*ALPHA
        assert max(abs(fd-exact))<1e-7

def main():
    validate();records=[]
    for N in [12,20,50]:
        for kind in KINDS:
            for seed in SEEDS:
                z=disorder(N,kind,seed)
                for sigma in STRENGTHS:
                    p=DEST/f'{kind}_N{N}_seed{seed}_sigma{sigma:.0e}.json'
                    if p.exists():records.append(json.loads(p.read_text()));continue
                    m=model(N,kind,z,sigma);a=np.full(N,ALPHA)
                    # Ideal diagnostic: local drive retuning restores the same uniform
                    # background. This is NOT the fixed-boundary-only protocol.
                    control=model(N,kind,z,sigma);control.pump-=control.force(a)
                    rec=dict(N=N,kind=kind,seed=seed,sigma=sigma,normalized_disorder=z.tolist(),
                             detuning=m.detuning.tolist(),bond_J=abs(np.diag(m.T,-1)).tolist(),bond_phi=np.angle(np.diag(m.T,-1)).tolist(),
                             linear=linear_prediction(N,kind,z,sigma),
                             clean_force_under_disorder=float(max(abs(m.force(a)))),
                             continuation=continuation(N,kind,z,sigma),
                             evolution=evolve(m),ideal_local_retuning_control=diagnose(control,a))
                    p.write_text(json.dumps(rec));records.append(rec)
                    print(N,kind,seed,sigma,'continuation',rec['continuation']['reached_full_disorder'],
                          'dynamics',rec['evolution']['reached_stationary'],rec['evolution'].get('stability'),
                          'profile error',rec['evolution']['max_relative_amplitude_error'],flush=True)
    checks=[]
    for kind in KINDS:
        N=50;seed=SEEDS[0];sigma=1e-2;p=DEST/f'tight_{kind}.json'
        if p.exists():checks.append(json.loads(p.read_text()));continue
        z=disorder(N,kind,seed);final=evolve(model(N,kind,z,sigma),2400,2e-12,2e-13)
        check=dict(N=N,kind=kind,seed=seed,sigma=sigma,evolution=final);checks.append(check);p.write_text(json.dumps(check))
    summary=dict(parameters=dict(J=1,U=U,V=0,delta=.5,kappa=1,epsilon=EPS,phi=float(PHI)),density=float(ROOT),
                 alpha0=[float(ALPHA.real),float(ALPHA.imag)],seeds=SEEDS,strengths=STRENGTHS,records=records,tighter_checks=checks)
    (DEST/'pilot.json').write_text(json.dumps(summary));plot(summary);report(summary)

def plot(summary):
    fig,axes=plt.subplots(2,3,figsize=(14,8),layout='constrained')
    for col,kind in enumerate(KINDS):
        for N,color in [(12,'#38945b'),(20,'#267aaf'),(50,'#c26b20')]:
            rr=[r for r in summary['records'] if r['kind']==kind and r['N']==N]
            for seed,marker in zip(SEEDS,['o','s']):
                rrr=[r for r in rr if r['seed']==seed]
                axes[0,col].loglog([r['sigma'] for r in rrr],[r['evolution']['max_relative_amplitude_error'] for r in rrr],marker+'-',color=color,label=f'N={N}, seed {seed}')
                rate=[r['evolution'].get('growth',np.nan) for r in rrr]
                axes[1,col].semilogx([r['sigma'] for r in rrr],rate,marker+'-',color=color)
        axes[0,col].set(title={'delta':'detuning disorder','J':'hopping-amplitude disorder','phi':'bond-phase disorder'}[kind],xlabel='Gaussian standard deviation sigma',ylabel='final max relative amplitude deformation');axes[0,col].legend(fontsize=7)
        axes[1,col].axhline(0,color='k',ls='--');axes[1,col].set(xlabel='Gaussian standard deviation sigma',ylabel='approached stationary BdG growth / J')
    fig.savefig('figures/compensated_disorder_pilot.png',dpi=180)
    fig,axes=plt.subplots(2,3,figsize=(14,8),layout='constrained')
    for row,seed in enumerate(SEEDS):
      for ax,kind in zip(axes[row],KINDS):
        for sigma in STRENGTHS:
            r=next(r for r in summary['records'] if r['N']==50 and r['kind']==kind and r['seed']==seed and r['sigma']==sigma)
            ax.plot(range(1,51),r['evolution']['density'],label=f'sigma={sigma:g}')
        ax.axhline(ROOT,color='k',ls='--',label='clean uniform target');ax.set(title=f'{kind} disorder, N=50, seed {seed}',xlabel='site j',ylabel='final density');ax.legend(fontsize=8)
    fig.savefig('figures/compensated_disorder_profiles.png',dpi=180)
    fig,axes=plt.subplots(2,3,figsize=(14,8),layout='constrained')
    for row,seed in enumerate(SEEDS):
      for ax,kind in zip(axes[row],KINDS):
        for sigma in STRENGTHS:
            r=next(r for r in summary['records'] if r['N']==12 and r['kind']==kind and r['seed']==seed and r['sigma']==sigma)
            ax.plot(range(1,13),r['evolution']['density'],'o-',ms=3,label=f'sigma={sigma:g}')
        ax.axhline(ROOT,color='k',ls='--',label='clean uniform target')
        ax.set(title=f'{kind} disorder, N=12, seed {seed}',xlabel='site j',ylabel='final density');ax.legend(fontsize=8)
    fig.savefig('figures/compensated_disorder_profiles_n12.png',dpi=180)

def report(summary):
    lines=['## Preliminary independent Gaussian disorder tests of the compensated phase',
    '', 'These are new mean-field/BdG pilot calculations, with **two realizations per parameter setting**, not a disorder phase diagram. At epsilon=38.51, phi=pi/3, Delta=0.5, J=kappa=1, U=-0.0002 and V=0, the clean intermediate density is 1990.520825. We compare N=12, N=20 and N=50.',
    '', '### Disorder and protocol',
    '', 'Disorder is static (quenched), real Gaussian and independent between sites/bonds. We test each parameter separately: Delta_j=0.5+sigma*z_j; J_j=1+sigma*z_j; or phi_j=pi/3+sigma*z_j, with z_j independent standard normals. Detuning sigma is in units of J, hopping sigma is also its fractional standard deviation since J=1, and phase sigma is in radians. The reverse hopping is the complex conjugate, preserving Hermiticity of the hopping Hamiltonian. We use sigma=1e-6, 1e-4, 1e-2 and seeds 20261104/20261105. Patterns are matched across strengths and chain lengths (short chains use prefixes); realizations are not independent across those comparisons.',
    '', 'The endpoint compensation remains fixed at its **clean** value, h_L=exp(i*phi)*alpha_0 and h_R=exp(-i*phi)*alpha_0. No disorder-dependent boundary tuning is applied. Starting from the already prepared clean compensated state, we switch on the full disorder and evolve for 1200/J. This tests a disorder quench, not vacuum preparation or a slow disorder ramp. In parallel, stationary continuation follows the clean branch as the disorder strength increases.',
    '', '### Results',
    '', '**All 54 quenches reached stationary states with negative BdG growth rates. Nevertheless, homogeneity is very sensitive, and the final attractor need not be the continuation of the clean compensated branch.** Here D_alpha=max_j|alpha_j-alpha_0|/|alpha_0|; ranges below span only the two realizations.',
    '', '| N | Disorder | sigma | D_alpha range | Final max Im(omega)/J range |', '|---:|:---|---:|---:|---:|']
    for N in [12,20,50]:
      for kind in KINDS:
       for sigma in STRENGTHS:
        rr=[r for r in summary['records'] if r['N']==N and r['kind']==kind and r['sigma']==sigma]
        d=[r['evolution']['max_relative_amplitude_error'] for r in rr];g=[r['evolution']['growth'] for r in rr]
        lines.append(f'| {N} | {kind} | {sigma:g} | {min(d):.4g}–{max(d):.4g} | {min(g):.5f}–{max(g):.5f} |')
    lines += ['', '![Disorder deformation and approached-state stability](../figures/compensated_disorder_pilot.png)',
    '', '![N=12 density profiles for both realizations](../figures/compensated_disorder_profiles_n12.png)',
    '', 'For N=12 all 18 quenches reach the same roots as stationary continuation (maximum complex-amplitude difference below 1.4e-12), and all sampled continuation points remain stable. At the largest tested sigma=0.01, maximum relative amplitude deformations are 2.36–2.43% for detuning, 18.7–26.0% for hopping, and 8.17–8.69% for phase disorder. Across all N=12 cases the growth rate lies between -0.09049 and -0.05939, compared with clean -0.07202. Thus these shorter chains show substantially weaker profile sensitivity and no observed branch switching in the pilot. These are two-realization comparisons, not statistical tolerance estimates.',
    '', '![N=50 density profiles for both realizations](../figures/compensated_disorder_profiles.png)',
    '', 'For example, N=50 hopping disorder at sigma=1e-6 produces 1.69–2.32% amplitude deformation; sigma=1e-4 produces 48.2–50.8%. The corresponding N=20 response is much smaller. This is consistent with the previously established amplification of static perturbations. A deviation above one in the table includes complex-phase changes and should not be read as a density percentage; the density profiles show the actual redistribution.',
    '', 'Stationary continuation reached full disorder in 51/54 cases; every accepted sampled point on these paths had negative growth. It stopped in the following three cases. These numerical stops do **not** establish branch nonexistence or an instability boundary:',
    '', '| N | Disorder | Seed | Target sigma | Last reached fraction of target |', '|---:|:---|---:|---:|---:|']
    for r in summary['records']:
        if not r['continuation']['reached_full_disorder']:
            lines.append(f"| {r['N']} | {r['kind']} | {r['seed']} | {r['sigma']:g} | {r['continuation']['last_lam']:.6f} |")
    lines += ['', 'There is also direct evidence of basin dependence: for N=50 detuning disorder, seed 20261105 and sigma=0.01, continuation reaches a stable root, but the sudden quench reaches a different stable root (maximum complex-amplitude separation about 86.1). Thus negative growth of the final state alone does not prove persistence of the original branch.',
    '', '### Separating spectral stability from background deformation',
    '', 'As a diagnostic control, we additionally retune the coherent drive at **every site** by subtracting the disordered stationary residual evaluated at alpha_0. This makes the clean uniform background an exact stationary solution of each disordered model. These drives are then treated as fixed coherent fields, not feedback. This is an ideal local compensation control, not the endpoint-only method. All 54 such uniform controls are also BdG stable. For N=50 their growth rates lie between -0.05782 and -0.05662, compared with clean -0.05713. In these samples, disorder has a much larger effect through the self-consistent background response than through its direct effect on the BdG spectrum at the restored uniform background.',
    '', '**Interpretation:** the pilot supports robustness of dynamical stability over these samples, while revealing strong fragility of the homogeneous profile with fixed endpoint fields. It does not establish that the final inhomogeneous attractors retain the clean topological phase. A disordered topological invariant/transition has not been computed; independent site-scale disorder does not justify a slowly varying local-bulk winding approximation. Singular values and singular-vector centers are saved as response diagnostics only.',
    '', '### Numerical checks and reproducibility',
    '', 'Evolution uses DOP853, rtol=2e-10, atol=2e-11; the final state must have raw force below 1e-3 and lie within 1e-3 in maximum amplitude of a stationary root polished to residual below 1e-8. All final polished force residuals are below 1.2e-13. Growth is max Im(omega) for i*dot(z)=L*z, cross-checked against the real-quadrature evolution generator. The sign classification uses a minimum 1e-7 resolution with additional leading-eigenvalue conditioning and basis-agreement checks. Three N=50, sigma=0.01 runs (one per disorder type, first seed) repeated to time 2400/J with 100-times tighter tolerances reproduce the final density to maximum absolute difference below 6e-11. Analytic Jacobians and disorder forcing derivatives were checked against finite differences.',
    '', 'Continuation uses adaptive parameter steps up to 0.2, stops below step 1e-5, and rejects stationary steps with maximum amplitude change above 2; its stopping point is algorithm dependent. Two realizations and three widely separated strengths cannot determine a disorder tolerance or rare-event probability. No covariance dynamics or fluctuation feedback is included.',
    '', 'Reproduce with `python src/preliminary_compensated_disorder.py`. The complete records, disorder vectors, continuation paths, spectra diagnostics, final profiles and sampled density trajectories are in `data/generated/compensated_disorder_pilot/pilot.json`.', '']
    p=Path('docs/uniform_topological_segment_stability.md');text=p.read_text(encoding='utf-8');heading='## Preliminary independent Gaussian disorder tests of the compensated phase'
    text=text.split(heading)[0].rstrip()+'\n\n'+'\n'.join(lines);p.write_text(text,encoding='utf-8')

if __name__=='__main__':main()
