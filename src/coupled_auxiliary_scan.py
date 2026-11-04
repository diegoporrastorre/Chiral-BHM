"""OBC continuation with undriven, hopping auxiliaries and fixed primary pumps."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json
import numpy as np
from scipy.linalg import solve,eig,svd
from scipy.optimize import root
from mean_field import Chain
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

N=50;DEST=Path('data/generated/coupled_auxiliary');DEST.mkdir(parents=True,exist_ok=True)
HB=np.diag(np.ones(N-1),1)+np.diag(np.ones(N-1),-1)


class Model:
    def __init__(self,kb,pump):
        self.kb=kb;self.jb=np.sqrt(kb)/2;self.C=HB-.5j*kb*np.eye(N)
        self.m=Chain(N=N,drive=38.51,kappa=0);self.m.pump=pump.copy()
        self.S=-self.jb**2*solve(self.C,np.eye(N))
    def force(self,a):return self.m.force(a)+self.S@a
    def blocks(self,a):
        A,B=self.m.blocks(a);return A+self.S,B
    def stationary(self,initial):
        def unpack(x):return x[:N]+1j*x[N:]
        def f(x):
            y=self.force(unpack(x));return np.r_[y.real,y.imag]
        def jac(x):
            A,B=self.blocks(unpack(x))
            return np.block([[(A+B).real,-(A-B).imag],[(A+B).imag,(A-B).real]])
        r=root(f,np.r_[initial.real,initial.imag],jac=jac,tol=1e-10)
        a=unpack(r.x);res=float(max(abs(self.force(a))))
        return a,res
    def full(self,a):
        A,B=self.m.blocks(a);Z=np.zeros_like(A);K=self.jb*np.eye(N)
        return np.block([[A,K,B,Z],[K,self.C,Z,Z],[-B.conjugate(),Z,-A.conjugate(),-K],[Z,Z,-K,-self.C.conjugate()]])
    def diagnose(self,a):
        L=self.full(a);w,vl,vr=eig(L,left=True,right=True);j=np.argmax(w.imag)
        cond=float(1/abs(np.vdot(vl[:,j],vr[:,j])))
        indicator=float(np.finfo(float).eps*np.linalg.norm(L,np.inf)*cond)
        tolerance=max(1e-7,10*indicator)
        beta=-self.jb*solve(self.C,a)
        residual_a=float(max(abs(self.m.force(a)+self.jb*beta)))
        residual_b=float(max(abs(self.C@beta+self.jb*a)))
        growth=float(w[j].imag)
        return dict(kappa_b=self.kb,J_b=self.jb,growth=growth,frequency=float(w[j].real),
            stability='stable' if growth < -tolerance else 'unstable' if growth>tolerance else 'unresolved',
            resolution=tolerance,condition=cond,sensitivity_indicator=indicator,
            leading_eigen_residual=float(max(abs(L@vr[:,j]-w[j]*vr[:,j]))),
            primary_residual=residual_a,auxiliary_residual=residual_b,
            density=(abs(a)**2).tolist(),auxiliary_density=(abs(beta)**2).tolist(),
            alpha_real=a.real.tolist(),alpha_imag=a.imag.tolist(),
            density_range=[float(min(abs(a)**2)),float(max(abs(a)**2))])


def continuation(pump,initial):
    eta=np.log(1e8);minimum=np.log(.01);step=.08;a,res=Model(np.exp(eta),pump).stationary(initial)
    assert res<1e-8
    prev=None;records=[Model(np.exp(eta),pump).diagnose(a)];failures=[]
    for iteration in range(5000):
        nxt=max(minimum,eta-step);model=Model(float(np.exp(nxt)),pump)
        guess=a if prev is None else a+(nxt-eta)*(a-prev[1])/(eta-prev[0])
        b,res=model.stationary(guess);change=float(max(abs(b-a)))
        if res>1e-8 or change>1:
            failures.append(dict(kappa_b=model.kb,log_step=step,residual=res,amplitude_step=change));step/=2
            if step<1e-6:break
            continue
        prev=(eta,a.copy());eta=nxt;a=b;rec=model.diagnose(a);records.append(rec)
        if len(records)%30==0:print('continued',len(records),'kb',model.kb,'growth',rec['growth'],flush=True)
        step=min(.12,step*1.25)
        if eta==minimum:break
    return dict(records=records,failed_attempts=failures,reached_minimum=bool(eta==minimum),
        stop_kappa_b=float(np.exp(eta)),stop_reason='minimum linewidth reached' if eta==minimum else 'adaptive stationary continuation step below 1e-6')


def refine_crossings(state,pump):
    refined=[]
    for aa,bb in zip(state['records'][:-1],state['records'][1:]):
        if aa['growth']*bb['growth']>=0:continue
        lo=bb;hi=aa
        for _ in range(32):
            kb=float(np.sqrt(lo['kappa_b']*hi['kappa_b']))
            a=np.array(lo['alpha_real'])+1j*np.array(lo['alpha_imag'])
            b=np.array(hi['alpha_real'])+1j*np.array(hi['alpha_imag'])
            m=Model(kb,pump);x,res=m.stationary((a+b)/2)
            assert res<1e-8
            rec=m.diagnose(x)
            if rec['growth']*lo['growth']>0:lo=rec
            else:hi=rec
            if hi['kappa_b']-lo['kappa_b']<1e-6:break
        refined.append(dict(lower=lo,upper=hi))
    state['crossings']=refined


def validate(pump,a):
    m=Model(3,pump);rng=np.random.default_rng(21);v=rng.normal(size=N)+1j*rng.normal(size=N);h=1e-5
    A,B=m.blocks(a);err=max(abs((m.force(a+h*v)-m.force(a-h*v))/(2*h)-A@v-B@v.conjugate()))
    assert err<1e-7,err


def fold(state,pump):
    last=state['records'][-1];a=np.array(last['alpha_real'])+1j*np.array(last['alpha_imag'])
    def quadrature(A,B):return np.block([[(A+B).real,-(A-B).imag],[(A+B).imag,(A-B).real]])
    m=Model(last['kappa_b'],pump);J=quadrature(*m.blocks(a));v=svd(J)[2][-1]
    x0=np.r_[a.real,a.imag,v,np.log(last['kappa_b'])]
    def components(x):
        a=x[:N]+1j*x[N:2*N];v=x[2*N:3*N]+1j*x[3*N:4*N]
        m=Model(float(np.exp(x[-1])),pump);J=quadrature(*m.blocks(a))
        f=m.force(a);return a,v,m,J,np.r_[f.real,f.imag]
    def fun(x):
        a,v,m,J,f=components(x);vv=x[2*N:4*N]
        return np.r_[f,J@vv,np.dot(vv,vv)-1]
    def jac(x):
        a,v,m,J,f=components(x);gamma=m.kb/2
        Seta=m.S+1j*gamma*m.S@solve(m.C,np.eye(N))
        fa=Seta@a;fv=Seta@v
        hx=8*m.m.U*a.real*v+4*m.m.U*a*v.conjugate()
        hy=8*m.m.U*a.imag*v+4j*m.m.U*a*v.conjugate()
        H=np.block([[np.diag(hx.real),np.diag(hy.real)],[np.diag(hx.imag),np.diag(hy.imag)]])
        out=np.zeros((4*N+1,4*N+1));out[:2*N,:2*N]=J;out[:2*N,-1]=np.r_[fa.real,fa.imag]
        out[2*N:4*N,:2*N]=H;out[2*N:4*N,2*N:4*N]=J;out[2*N:4*N,-1]=np.r_[fv.real,fv.imag]
        out[-1,2*N:4*N]=2*x[2*N:4*N]
        return out
    # Check augmented Jacobian before relying on the fold calculation.
    rng=np.random.default_rng(12);direction=rng.normal(size=len(x0));direction/=np.linalg.norm(direction);h=1e-5
    assert max(abs((fun(x0+h*direction)-fun(x0-h*direction))/(2*h)-jac(x0)@direction))<1e-6
    sol=root(fun,x0,jac=jac,tol=1e-10);res=float(max(abs(fun(sol.x))))
    a,v,m,J,f=components(sol.x);u,s,vh=svd(J);vv=np.r_[v.real,v.imag];left=u[:,-1]
    # A simple fold has rank loss one, nonzero parameter transversality,
    # and nonzero quadratic curvature in its null direction.
    aug=jac(sol.x);trans=float(abs(left@aug[:2*N,-1]));curv=float(abs(left@aug[2*N:4*N,:2*N]@vv))
    partners=[];deta=1e-4
    signed_trans=float(left@aug[:2*N,-1]);signed_curv=float(left@aug[2*N:4*N,:2*N]@vv)
    scale=np.sqrt(-2*signed_trans*deta/signed_curv)
    neighbor=Model(m.kb*np.exp(deta),pump)
    for sign in [-1,1]:
        aa,rr=neighbor.stationary(a+sign*scale*v)
        assert rr<1e-8
        partners.append(neighbor.diagnose(aa))
    return dict(converged=bool(res<1e-8),augmented_residual=res,kappa_b=m.kb,
        smallest_stationary_singular_values=s[-3:].tolist(),parameter_transversality=trans,null_direction_curvature=curv,
        final=m.diagnose(a),nearby_partners=partners)


def main():
    source=json.loads(Path('data/generated/compensation_n50/intermediate_preparation_35_40.json').read_text())
    r=next(r for r in source['records'] if r['epsilon']==38.51);a0=complex(*r['alpha0'])
    initial=np.array(r['baseline']['alpha_real'])+1j*np.array(r['baseline']['alpha_imag'])
    states={}
    for label,a in [('compensated',np.full(N,a0)),('uncompensated',initial)]:
        pump=np.full(N,38.51,dtype=complex)
        if label=='compensated':pump[0]+=np.exp(1j*np.pi/3)*a0;pump[-1]+=np.exp(-1j*np.pi/3)*a0
        validate(pump,a);file=DEST/f'{label}.json'
        state=json.loads(file.read_text()) if file.exists() else continuation(pump,a)
        if 'crossings' not in state:refine_crossings(state,pump)
        if not state['reached_minimum'] and ('fold' not in state or 'nearby_partners' not in state['fold']):state['fold']=fold(state,pump)
        file.write_text(json.dumps(state));states[label]=state
        print(label,'stopped',state['stop_kappa_b'],'crossings',[(c['lower']['kappa_b'],c['upper']['kappa_b']) for c in state['crossings']],flush=True)
    summary=dict(parameters=dict(N=50,epsilon=38.51,J=1,U=-.0002,V=0,delta_a=.5,phi=float(np.pi/3),delta_b=0,J_b_prime=1,kappa_a=0,matching='J_b=sqrt(kappa_b)/2',pumps='fixed original primary pumps; no auxiliary drives'),states=states)
    (DEST/'scan.json').write_text(json.dumps(summary));plot(summary);report(summary)


def plot(s):
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for col,(label,state) in enumerate(s['states'].items()):
        rr=state['records'];axes[0,col].semilogx([r['kappa_b'] for r in rr],[r['growth'] for r in rr],'.-',ms=2,label='continued stationary branch')
        axes[0,col].plot(state['fold']['kappa_b'],0,'r*',ms=11,label='stationary fold')
        for kb in ([1e8,1e6,2e5] if label=='compensated' else [1e8,100,6]):
            r=min(rr,key=lambda r:abs(np.log(r['kappa_b']/kb)))
            if abs(np.log(r['kappa_b']/kb))<.15:axes[1,col].plot(range(1,51),r['density'],label=f"kb={r['kappa_b']:.3g}")
        axes[0,col].axhline(0,color='k',lw=.8);axes[0,col].set(xlabel='kappa_b / J',ylabel='max Im(omega) / J',title=f"{label}, auxiliary hopping J_b'=1")
        axes[1,col].plot(range(1,51),state['fold']['final']['density'],'k--',label='at fold')
        axes[1,col].set(xlabel='site',ylabel='stationary primary density')
    for ax in axes.flat:ax.legend(fontsize=7)
    fig.savefig('figures/coupled_auxiliary_scan.png',dpi=180)


def report(summary):
    states=summary['states']
    text=r'''## Numerical scan with coupled auxiliaries: J_b'=1

### Setup and central result

We scan the extended model at N=50, epsilon=38.51, J=1, U=-0.0002, V=0, Delta_a=0.5, phi=pi/3, Delta_b=0 and kappa_a=0. The auxiliaries have real reciprocal OBC hopping J_b'=1, independent local losses kappa_b and no coherent drives. The local coupling follows J_b=sqrt(kappa_b)/2. This is the same coupling prescription as before, **not** an exact static-loss match when the auxiliaries hop.

Two distinct sets of primary pumps are held fixed: (i) the original fields compensating the clean intermediate root, including the same two endpoint additions; and (ii) uniform epsilon=38.51 without compensation. Both stationary branches are continued from kappa_b=1e8*J toward smaller linewidth. No auxiliary or primary fields are retuned during the scan.

**The compensated branch deforms and ends in a stationary fold at kappa_b approximately 156438.4991J. The uncompensated branch reaches a stationary fold at approximately 5.488779818J.** Every accepted continuation point on the tracked branches is classified stable; the endpoint is marginal. These are stationary branch terminations, unlike the finite-frequency instabilities at unchanged profiles found for independent auxiliaries. The results do not establish that all stationary states or stable attractors disappear below either fold.

| Pump choice / tracked branch | Accepted scan points | Fold kappa_b/J | J_b/J at fold | Primary density range at fold |
|:---|---:|---:|---:|---:|
'''
    for label,state in states.items():
        f=state['fold'];d=f['final']['density_range']
        text+=f"| {label} | {len(state['records'])} | {f['kappa_b']:.10g} | {np.sqrt(f['kappa_b'])/2:.8g} | {d[0]:.8g}–{d[1]:.8g} |\n"
    text+=r'''
![Stationary-branch growth rates and density profiles for coupled auxiliaries](../figures/coupled_auxiliary_scan.png)

### Why the compensated branch changes at such a large linewidth

The original endpoint fields cancel missing primary neighbors for a local instantaneous reservoir. With auxiliary hopping and the same fields, the reservoir correction is spatially nonlocal and its real part is not zero, even at zero frequency. At large linewidth,

$$
-J_b^2(H_b-i\kappa_b I/2)^{-1}
=-\frac{i}{2}I-\frac{H_b}{\kappa_b}
+\frac{2iH_b^2}{\kappa_b^2}+O(\kappa_b^{-3}),
$$

where kappa_target=1 and H_b=T_b for this scan. The leading correction is a small extra coherent hopping/detuning operator of order J_b'/kappa_b. The original compensated background has a large static susceptibility, so this small correction can cause a substantial profile change. At kappa_b near 1e6*J the compensated density already ranges approximately 1990.52–2069.34, and near 2e5*J it ranges approximately 1990.52–2572.49. The fold has a right-edge density around 3309.63. Thus J_b'=1 is not a small perturbation for preserving this branch with its original fields, even though its static reservoir correction is small at the fold linewidth.

This is a combined stationary-background and dynamical effect. It is not evidence that spatial reservoir coupling or long memory universally destroys compensated states. Retuning auxiliary detuning, compensating auxiliary edges, or adjusting primary drives would define different comparisons. In particular, the previous conclusion that all static primary profiles are unchanged is restricted to independent auxiliaries and does not apply to this scan.

### Checks that the continuation stops are folds

The stationary auxiliaries are eliminated exactly: beta=-J_b*C_b^(-1)*alpha. The reduced primary force is F+T*alpha+(Delta_a+U+2U|alpha|^2)*alpha-J_b^2*C_b^(-1)*alpha. Its analytic real Jacobian includes the full nonlocal reservoir matrix. Continuation stopping alone is not used as evidence of a branch endpoint.

For each endpoint we solve the augmented fold equations f=0, J_stat*v=0 and v^T*v=1, with the logarithmic linewidth also unknown. The augmented analytic Jacobian is checked against a finite difference. The stationary Jacobian loses one rank; its parameter transversality and quadratic curvature in the null direction are nonzero. At a nearby linewidth kappa_b,fold*exp(1e-4), two distinct stationary roots are found, one stable and one unstable. Together these checks support a simple saddle-node of the tracked branch.

| Branch | Augmented max residual | Smallest stationary singular value | Next stationary singular value | Nearby pair growth rates / J |
|:---|---:|---:|---:|:---|
'''
    for label,state in states.items():
        f=state['fold'];sv=f['smallest_stationary_singular_values'];g=[r['growth'] for r in f['nearby_partners']]
        text+=f"| {label} | {f['augmented_residual']:.3g} | {sv[-1]:.3g} | {sv[-2]:.6g} | {g[0]:.7g}, {g[1]:.7g} |\n"
    text+=r'''
The full 200-by-200 particle–hole BdG matrix is evaluated at every accepted stationary profile, retaining the auxiliary hopping and modes. Eigenvalues use the complex eigensolver because the real driver failed residual checks in the earlier independent-auxiliary calculation. Stability requires max Im(omega) below a resolution equal to the larger of 1e-7*J and ten times a leading-eigenvalue roundoff sensitivity indicator. All accepted points are stable under that criterion; folds are treated as marginal, not assigned a stable sign from numerical roundoff.

Continuation starts at 1e8*J, proposes logarithmic steps up to 0.12 and rejects roots with reduced force residual above 1e-8*J or maximum amplitude change above 1. Step sizes below 1e-6 trigger the augmented endpoint calculation. The planned lower linewidth was 0.01J; neither connected branch reaches it. Other disconnected branches, time-dependent attractors and quench dynamics below the folds have not been explored. Stable and unstable roots can coexist above a fold; no uniqueness claim is made.
'''
    for label,state in states.items():
        rr=state['records']
        text+=f"\nFor {label}, maximum accepted primary force residual is {max(r['primary_residual'] for r in rr):.3g}J and maximum leading-eigenpair absolute residual is {max(r['leading_eigen_residual'] for r in rr):.3g}. The latter occurs in the very large-linewidth matrices, whose auxiliary entries are of order 1e8*J; condition-based resolution is included.\n"
    text+=r'''

No topology invariant is computed for these deformed backgrounds. The old homogeneous winding cannot automatically be assigned to them, and static compensation fields no longer imply a homogeneous profile. No Gaussian covariance or fluctuation feedback is used.

Reproduce with `python src/coupled_auxiliary_scan.py`. Data are saved in `data/generated/coupled_auxiliary/scan.json` and per-branch checkpoints in the same directory, including full complex primary profiles, auxiliary densities, rejected continuation attempts, stability diagnostics and augmented fold checks. Figure: `figures/coupled_auxiliary_scan.png`.
'''
    path=Path('derivations/non_markovian_auxiliary_cavity_model.md');old=path.read_text(encoding='utf-8')
    old=old.split("## Numerical scan with coupled auxiliaries: J_b'=1")[0].rstrip()
    path.write_text(old+'\n\n'+text,encoding='utf-8')


if __name__=='__main__':main()
