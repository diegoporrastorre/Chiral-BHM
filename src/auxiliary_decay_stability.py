"""Matched resonant auxiliary-decay stability, fixed V=0 mean-field profiles."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json
import numpy as np
from scipy.linalg import eig,eigvals
from scipy.optimize import brentq
from mean_field import Chain
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DEST=Path('data/generated/auxiliary_decay');DEST.mkdir(parents=True,exist_ok=True)
N=50;TARGET=1.


def matrices(a):
    m=Chain(N=N,drive=38.51,kappa=0)
    A,B=m.blocks(a);L=np.block([[A,B],[-B.conjugate(),-A.conjugate()]])
    return m,A,B,L


def full(A,B,kb):
    I=np.eye(N);Z=np.zeros_like(A);jb=np.sqrt(TARGET*kb)/2;C=-.5j*kb*I
    return np.block([[A,jb*I,B,Z],[jb*I,C,Z,Z],[-B.conjugate(),Z,-A.conjugate(),-jb*I],[Z,Z,-jb*I,-C.conjugate()]])


def modes(lam,kb):
    gamma=kb/2;jb2=TARGET*kb/4
    # Stable quadratic evaluation: compute large root, obtain small root by product.
    discr=np.sqrt((lam+1j*gamma)**2+4*jb2)
    a=(lam-1j*gamma+discr)/2;b=(lam-1j*gamma-discr)/2
    product=-1j*gamma*lam-jb2
    mask=abs(a)>=abs(b)
    big=np.where(mask,a,b);small=product/big
    return np.r_[big,small]


def analyze(a,compensated):
    m,A,B,L=matrices(a);lam,left,right=eig(L,left=True,right=True)
    conditions=1/abs(np.sum(left.conjugate()*right,axis=0))
    y=lam.imag;x=lam.real;h=TARGET/2
    assert max(y)<h
    thresholds=np.zeros(len(lam));positive=y>0
    thresholds[positive]=2*(y[positive]+abs(x[positive])*np.sqrt(y[positive]/(h-y[positive])))
    idx=int(np.argmax(thresholds));crit=float(thresholds[idx])
    grid=np.unique(np.r_[np.logspace(-4,4,401),crit*np.array([.99,1,1.01]),[.01,.1,1,2,5,10,20,100,1000]])
    records=[]
    for kb in grid:
        w=modes(lam,kb);j=int(np.argmax(w.imag))
        records.append(dict(kappa_b=float(kb),J_b=float(np.sqrt(kb)/2),growth=float(w[j].imag),frequency=float(w[j].real),unstable_modes=int(sum(w.imag>1e-8))))
    def growth(kb):return float(max(modes(lam,kb).imag))
    refined=brentq(growth,crit*.99,crit*1.01,xtol=1e-11)
    checks=[]
    # Full 4N diagonalization and real-quadrature generator independently check mapping.
    for kb in [.01,.1,1,crit*.99,crit,crit*1.01,10,100,1000]:
        Lfull=full(A,B,kb);w,vl,vr=eig(Lfull,left=True,right=True);j=np.argmax(w.imag)
        predicted=modes(lam,kb)
        distance=float(max(min(abs(v-predicted)) for v in w))
        AA=np.block([[A,np.sqrt(kb)/2*np.eye(N)],[np.sqrt(kb)/2*np.eye(N),-.5j*kb*np.eye(N)]])
        BB=np.block([[B,np.zeros_like(B)],[np.zeros_like(B),np.zeros_like(B)]])
        G=np.block([[(AA+BB).imag,(AA-BB).real],[-(AA+BB).real,(AA-BB).imag]])
        # Use complex LAPACK for the real-coordinate matrix too. The installed
        # real nonsymmetric driver returned spurious modes with large residuals
        # on these 200x200 matrices; casting preserves G while avoiding that path.
        wg,vg=eig(G.astype(complex));realgrowth=float(max(wg.real))
        real_residual=float(np.max(abs(G@vg-vg*wg)))
        full_residual=float(np.max(abs(Lfull@vr-vr*w)))
        assert max(real_residual,full_residual)<1e-9
        cond=float(1/abs(np.vdot(vl[:,j],vr[:,j])))
        sensitivity=float(np.finfo(float).eps*np.linalg.norm(Lfull,np.inf)*cond)
        beta=-1j*a/np.sqrt(kb)
        pumps=np.full(N,38.51,dtype=complex)
        if compensated:pumps[0]+=np.exp(1j*np.pi/3)*a[0];pumps[-1]+=np.exp(-1j*np.pi/3)*a[0]
        primary=pumps+m.T@a+(.4998-.0004*abs(a)**2)*a+np.sqrt(kb)/2*beta
        aux=-.5j*kb*beta+np.sqrt(kb)/2*a
        checks.append(dict(kappa_b=float(kb),full_growth=float(w[j].imag),mapped_growth=growth(kb),real_growth=realgrowth,
            max_spectrum_mapping_distance=distance,leading_condition=cond,sensitivity_indicator=sensitivity,
            real_coordinate_eigen_residual=real_residual,full_eigen_residual=full_residual,
            primary_residual=float(max(abs(primary))),auxiliary_residual=float(max(abs(aux)))))
    return dict(critical_kappa_b=crit,critical_kappa_b_brent=refined,critical_J_b=float(np.sqrt(crit)/2),critical_memory_time=2/crit,
        triggering_lambda=[float(lam[idx].real),float(lam[idx].imag)],
        markov_growth=float(max(y)-.5),primary_max_eigen_condition=float(max(conditions)),
        primary_eigenvalues_real=lam.real.tolist(),primary_eigenvalues_imag=lam.imag.tolist(),
        density=(abs(a)**2).tolist(),records=records,checks=checks)


def main():
    source=json.loads(Path('data/generated/compensation_n50/intermediate_preparation_35_40.json').read_text())
    r=next(r for r in source['records'] if r['epsilon']==38.51)
    a0=complex(*r['alpha0']);a=np.array(r['baseline']['alpha_real'])+1j*np.array(r['baseline']['alpha_imag'])
    states={}
    for label,alpha,comp in [('compensated',np.full(N,a0),True),('uncompensated',a,False)]:
        states[label]=analyze(alpha,comp);s=states[label]
        print(label,'critical kb',s['critical_kappa_b'],'Jb',s['critical_J_b'],'lambda',s['triggering_lambda'],
            'max basis difference',max(abs(c['full_growth']-c['real_growth']) for c in s['checks']),flush=True)
    result=dict(parameters=dict(N=50,epsilon=38.51,J=1,U=-.0002,V=0,delta=.5,phi=float(np.pi/3),delta_b=0,kappa_a=0,kappa_target=1),states=states)
    (DEST/'stability.json').write_text(json.dumps(result,indent=2));plot(result);report(result)


def plot(s):
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    for label,state in s['states'].items():
        rr=state['records'];kb=np.array([r['kappa_b'] for r in rr]);g=np.array([r['growth'] for r in rr])
        axes[0].semilogx(kb,g,label=label);axes[0].axhline(state['markov_growth'],ls=':',alpha=.6)
        crit=state['critical_kappa_b'];mask=(kb>crit*.5)&(kb<crit*2)
        axes[1].plot(kb[mask],g[mask],'.-',label=label);axes[1].axvline(crit,ls='--',alpha=.5)
    for ax in axes:ax.axhline(0,color='k',lw=.8);ax.legend();ax.set(xlabel='auxiliary linewidth kappa_b / J',ylabel='max Im(omega) / J')
    axes[0].set(title='Matched static loss: 4 J_b^2 / kappa_b = 1');axes[1].set(title='Stability crossings')
    fig.savefig('figures/auxiliary_decay_stability.png',dpi=180)


def report(result):
    states=result['states']
    lines=[r'''# Stability with matched non-Markovian auxiliary decay

## Result

For N=50, epsilon=38.51, J=1, U=-0.0002, V=0, Delta_a=0.5 and phi=pi/3, replacing direct kappa=1 loss by resonant local auxiliaries can destabilize both previously stable stationary backgrounds. We use Delta_b=0, no direct primary loss, and J_b=sqrt(kappa_b)/2, so 4*J_b^2/kappa_b=1. All primary stationary profiles remain exactly the same as in the Markovian problem; only dynamics and auxiliary populations change.

The compensated homogeneous intermediate state is stable for kappa_b>7.41212493*J. The uncompensated stationary state previously reached from vacuum is stable for kappa_b>5.45219985*J. Equality is marginal, and smaller positive linewidths are unstable for these respective branches. These are local linear-stability statements about fixed stationary states, not results of evolving from vacuum in the new model.

| Background | Critical kappa_b/J | J_b/J at crossing | Critical memory time J*tau_b | Markovian growth/J |
|:---|---:|---:|---:|---:|''']
    for label,s in states.items():
        lines.append(f"| {label} | {s['critical_kappa_b']:.9f} | {s['critical_J_b']:.9f} | {s['critical_memory_time']:.9f} | {s['markov_growth']:.9f} |")
    lines += [r'''
Here tau_b=2/kappa_b. Thus increasing memory from the Markovian limit first destabilizes the compensated state. For 5.4522<kappa_b/J<7.4121 the uncompensated branch remains stable while the compensated branch is unstable; below 5.4522 both tested branches are unstable. This does not exclude additional stable stationary branches or nonstationary attractors.

![Auxiliary linewidth scan and stability crossings](../figures/auxiliary_decay_stability.png)

## Numerical growth rates

Growth means max Im(omega) for i*dot(z)=L_ab*z. A negative value is asymptotic linear stability.

| kappa_b/J | J_b/J | Compensated growth/J | Uncompensated growth/J |
|---:|---:|---:|---:|''']
    for kb in [.01,.1,1,2,5,10,20,100,1000]:
        g=[next(r['growth'] for r in states[label]['records'] if r['kappa_b']==kb) for label in ['compensated','uncompensated']]
        lines.append(f'| {kb:g} | {np.sqrt(kb)/2:.6g} | {g[0]:.9f} | {g[1]:.9f} |')
    lines += [r'''
The large-linewidth limit recovers the old Markovian growth rates. Convergence need not be monotonic: finite memory can slightly increase damping before longer memory causes instability. Decreasing kappa_b at fixed J_b would be a different comparison because it changes the stationary effective loss.

## Exact spectral mapping and threshold

Let L_a^(0) denote the primary-chain 2N-by-2N BdG operator at the chosen background, with the direct loss removed. Its blocks retain the onsite U*n^2 convention, including +U. Let lambda=x+i*y be any eigenvalue of L_a^(0). Since the resonant auxiliary particle and hole blocks are both -i*gamma_b times the identity, and the squared primary–auxiliary Nambu coupling is J_b^2 times the identity, the enlarged characteristic polynomial factorizes into scalar quadratic factors:

$$
(\omega-\lambda)(\omega+i\gamma_b)-J_b^2=0,
\qquad \gamma_b=\kappa_b/2,\quad J_b^2=\kappa_{\rm target}\gamma_b/2.
$$

Each lambda therefore gives two eigenvalues

$$
\omega_\pm=\frac{\lambda-i\gamma_b
\pm\sqrt{(\lambda+i\gamma_b)^2+4J_b^2}}{2}.
$$

This relation accounts for all 4N modes, including the auxiliary-like modes. The original Markovian spectrum is simply lambda-i*kappa_target/2. Both backgrounds have max y<kappa_target/2 but include modes with y>0: removing instantaneous damping exposes their gain, which the finite-bandwidth auxiliary reservoir may fail to suppress.

For y<=0, coupling this scalar mode to a damped auxiliary does not generate an instability. For 0<y<h, h=kappa_target/2, the two roots are strictly damped precisely when

$$
\gamma_b>y+|x|\sqrt{\frac{y}{h-y}}.
$$

Thus the full-chain threshold is

$$
\boxed{\kappa_{b,c}=2\max_{\lambda:\,0<y<h}
\left[y+|x|\sqrt{\frac{y}{h-y}}\right].}
$$

For nonzero x, this follows by setting a marginal omega real: the imaginary part of the quadratic gives omega=gamma_b*x/(gamma_b-y), and the real part gives y[1+x^2/(gamma_b-y)^2]=h. Stability requires gamma_b>y, selecting the upper crossing. The x=0 limit gives gamma_b,c=y. The high-linewidth branch is stable, and each scalar pair has only this stability boundary on gamma_b>y. Hence the matched family has no hidden reentrant stable interval below the reported threshold for these fixed backgrounds.

The modes that set the thresholds are:

| Background | Triggering primary eigenvalue lambda/J | Marginal frequency magnitude / J |
|:---|:---|---:|''']
    for label,s in states.items():
        x,y=s['triggering_lambda'];gamma=s['critical_kappa_b']/2
        freq=abs(gamma*x/(gamma-y))
        lines.append(f'| {label} | {x:.9f} + {y:.9f} i | {freq:.9f} |')
    lines += [r'''
Particle–hole symmetry supplies the opposite-frequency partner. The loss of stability is oscillatory (a finite-frequency linear crossing), not a stationary root disappearing. A nonlinear limit cycle is a possible subsequent behavior, but it has not been demonstrated here. The triggering mode need not be the least-damped mode in the Markovian limit; its nonzero frequency and the frequency-dependent reservoir response matter.

## Topology and stationary profiles

The compensated density remains n=1990.52082476 on every primary site. Its original endpoint compensation is held fixed, and beta_j=-i*sqrt(kappa_target/kappa_b)*alpha_j. The uncompensated primary profile is the previously computed vacuum attractor at epsilon=38.51, with no endpoint additions; it is inhomogeneous and already contains an extended nonzero frozen-winding segment. We are not comparing a trivial bulk with a topological bulk in this study.

At omega_ref=0, the matched primary Schur complement is exactly the original Markovian BdG operator for either profile. The homogeneous compensated bulk therefore retains its point-gap winding -1 along the entire matched family even after its dynamics becomes unstable. No point-gap closure at zero is required for this finite-frequency instability. **Preservation of zero-frequency topology and existence of a stationary profile do not guarantee dynamical stability.** The winding statement is a bulk statement for the homogeneous target; no quantized invariant is assigned to the inhomogeneous chain here.

## Methods and checks

The logarithmic scan contains 401 base linewidth values from 1e-4*J to 1e4*J, plus representative points and values bracketing each threshold by 1%. The analytic threshold is checked with a Brent root of the growth function (absolute linewidth tolerance 1e-11*J). The quadratic is evaluated using the root product to avoid cancellation when auxiliary and primary time scales are widely separated. Nine representative points per background, including both sides of each crossing, are independently diagonalized as full 200-by-200 matrices and in real-quadrature coordinates.

The installed real nonsymmetric eigensolver returned spurious eigenpairs with large residuals for these 200-dimensional generators. Those results are rejected. The real-coordinate matrices are instead cast to complex without changing any entries and diagonalized with the complex eigensolver; every accepted eigenpair set is checked against the original matrix. This agrees with the independent quadratic spectral mapping. This numerical issue is not interpreted as a physical instability.
''']
    for label,s in states.items():
        c=s['checks']
        lines.append(f"For {label}: maximum full-spectrum mapping distance {max(r['max_spectrum_mapping_distance'] for r in c):.3g}J; maximum growth disagreement between particle/hole and real coordinates {max(abs(r['full_growth']-r['real_growth']) for r in c):.3g}J; maximum full eigenpair residual {max(r['full_eigen_residual'] for r in c):.3g}; primary/auxiliary stationary force residuals below {max(max(r['primary_residual'],r['auxiliary_residual']) for r in c):.3g}J.")
    lines += [r'''
The sign classification away from crossings uses a descriptive 1e-7*J resolution; the 1% brackets are well outside this band. Leading-eigenvalue condition numbers and first-order roundoff indicators are saved. The reported digits describe this N=50 finite-dimensional mean-field/BdG model, not a thermodynamic boundary or an interval-arithmetic certificate. No time evolution, basin study or fluctuation/noise calculation is performed, and no covariance feeds back into the stationary background.

Model and analytical conventions: [auxiliary-cavity derivation](../derivations/non_markovian_auxiliary_cavity_model.md). Reproduce with `python src/auxiliary_decay_stability.py`. Results: `data/generated/auxiliary_decay/stability.json`; figure: `figures/auxiliary_decay_stability.png`. Primary backgrounds are taken from the epsilon=38.51 entry in `data/generated/compensation_n50/intermediate_preparation_35_40.json`.
''']
    Path('docs/auxiliary_decay_stability.md').write_text('\n\n'.join(lines),encoding='utf-8')


if __name__=='__main__':main()
