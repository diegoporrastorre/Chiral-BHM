"""Local-detuning homodyne sensing at fixed mean-field backgrounds (no feedback)."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import json
from pathlib import Path
import numpy as np
from scipy.linalg import solve, eigvals, solve_continuous_lyapunov
from mean_field import Chain
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DEST=Path('data/generated/compensated_sensing');DEST.mkdir(parents=True,exist_ok=True)
N=50;KAPPA=1.


def build(compensated,a0):
    m=Chain(N=N,drive=32);m.pump=m.pump.astype(complex)
    if compensated:
        m.pump[0]+=np.exp(1j*m.phi)*a0;m.pump[-1]+=np.exp(-1j*m.phi)*a0
    return m


def analyze(m,a):
    A,B=m.blocks(a)
    G=np.block([[(A+B).imag,(A-B).real],[-(A+B).real,(A-B).imag]])
    F=np.vstack([np.diag(np.sqrt(2)*a.imag),np.diag(-np.sqrt(2)*a.real)])
    R=-solve(G,F);D=-np.sqrt(KAPPA)*R
    M=np.eye(2*N)+KAPPA*solve(G,np.eye(2*N))
    C=.5*M@M.T
    V=solve_continuous_lyapunov(G,-.5*KAPPA*np.eye(2*N))
    occ=(np.diag(V)[:N]+np.diag(V)[N:]-1)/2
    results=[]
    for ell in [0,N-1]:
        ix=[ell,N+ell];c=C[np.ix_(ix,ix)]
        for j in range(N):
            d=D[ix,j];u=solve(c,d);u/=np.linalg.norm(u)
            if u@d<0:u=-u
            ug=d/np.linalg.norm(d)
            noise=float(u@c@u);slope=float(u@d)
            results.append(dict(signal_site=j+1,output_site=ell+1,theta=float(np.arctan2(u[1],u[0])),
                u=u.tolist(),slope=slope,noise_psd=noise,detuning_floor=np.sqrt(noise)/slope,
                max_signal_slope=float(np.linalg.norm(d)),
                gain_quadrature_noise=float(ug@c@ug),gain_quadrature_floor=float(np.sqrt(ug@c@ug)/np.linalg.norm(d)),
                output_covariance=c.tolist()))
    omega=np.block([[np.zeros((N,N)),np.eye(N)],[-np.eye(N),np.zeros((N,N))]])
    return dict(growth=float(eigvals(G).real.max()),residual=float(max(abs(m.force(a)))),
        population=float(sum(abs(a)**2)),edge_flux=(KAPPA*abs(a[[0,-1]])**2).tolist(),
        total_coherent_emission=float(KAPPA*sum(abs(a)**2)),
        fluctuation_occupancy=occ.tolist(),max_fluctuation_to_coherent=float(max(occ/abs(a)**2)),
        lyapunov_relative_residual=float(np.linalg.norm(G@V+V@G.T+.5*KAPPA*np.eye(2*N))/max(1,np.linalg.norm(V))),
        symplectic_absolute_error=float(np.max(abs(M@omega@M.T-omega))),
        best=min(results,key=lambda r:r['detuning_floor']),measurements=results,
        response_real=(R[:N]/np.sqrt(2)).tolist(),response_imag=(R[N:]/np.sqrt(2)).tolist())


def nonlinear(m,a,analysis,j,output_site):
    rec=next(r for r in analysis['measurements'] if r['signal_site']==j+1 and r['output_site']==output_site)
    ell=rec['output_site']-1;u=np.array(rec['u'])
    deriv=np.array(analysis['response_real'])[:,j]+1j*np.array(analysis['response_imag'])[:,j]
    rows=[]
    for sign in [-1,1]:
        for amp in np.logspace(-10,-1,19):
            delta=float(sign*amp);mm=build(False,0)
            # Preserve precisely the original pumps and all model coefficients.
            mm.pump=m.pump.copy();mm.detuning=m.detuning.copy();mm.detuning[j]+=delta
            rr=mm.stationary(a+delta*deriv,1e-9);x=rr['alpha']
            actual=-np.sqrt(2*KAPPA)*np.array([(x[ell]-a[ell]).real,(x[ell]-a[ell]).imag])@u
            predicted=rec['slope']*delta
            growth=float(eigvals(mm.bogoliubov(x)).imag.max())
            rows.append(dict(delta=delta,converged=rr['converged'],residual=rr['residual'],growth=growth,
                relative_profile_change=float(max(abs(x-a))/max(abs(a))),
                readout=actual,prediction=predicted,relative_readout_error=float(abs(actual-predicted)/abs(predicted)),
                max_relative_amplitude_prediction_error=float(max(abs(x-a-delta*deriv))/max(abs(delta*deriv)))))
    return dict(signal_site=j+1,output_site=ell+1,theta=rec['theta'],records=rows)


def validate():
    # Passive vacuum scatterer: output remains vacuum, including direct input interference.
    m=Chain(N=N,U=0,drive=0);z=analyze(m,np.ones(N,dtype=complex))
    for r in z['measurements']:np.testing.assert_allclose(r['output_covariance'],.5*np.eye(2),atol=1e-12)
    assert max(abs(np.array(z['fluctuation_occupancy'])))<1e-12


def main():
    validate()
    src=json.loads(Path('data/generated/compensation_n50/preparation_30_35/epsilon_32.000000.json').read_text())
    a0=complex(*src['alpha0']);states={}
    for label,comp in [('uncompensated',False),('compensated',True)]:
        m=build(comp,a0)
        initial=np.full(N,a0) if comp else np.array(src['baseline']['alpha_real'])+1j*np.array(src['baseline']['alpha_imag'])
        rr=m.stationary(initial,1e-9);assert rr['converged'];a=rr['alpha']
        result=analyze(m,a);bestj=result['best']['signal_site']-1
        result['nonlinear']=[nonlinear(m,a,result,j,ell) for j in sorted({bestj,24}) for ell in [1,50]]
        # Independent complex particle/hole static response check.
        L=m.bogoliubov(a);forcing=np.vstack([np.diag(a),-np.diag(a.conjugate())])
        complex_response=-solve(L,forcing)[:N]
        real_response=np.array(result['response_real'])+1j*np.array(result['response_imag'])
        result['response_basis_relative_error']=float(np.max(abs(complex_response-real_response))/np.max(abs(real_response)))
        assert result['response_basis_relative_error']<1e-9
        # Tighter stationary comparison is recorded by the nonlinear force residuals.
        states[label]=result
        print(label,'best',result['best'],'population',result['population'],'max fluct ratio',result['max_fluctuation_to_coherent'],flush=True)
    summary=dict(N=N,epsilon=32,phi=float(np.pi/3),delta=.5,U=-.0002,V=0,kappa=1,J=1,
        convention='two-sided symmetrized PSD, vacuum quadrature PSD=1/2',states=states)
    (DEST/'sensing.json').write_text(json.dumps(summary,indent=2));plot(summary);report(summary)


def plot(s):
    fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    for label,state in s['states'].items():
        for ell,ls in [(1,'--'),(50,'-')]:
            rr=[r for r in state['measurements'] if r['output_site']==ell]
            axes[0,0].semilogy(range(1,51),[r['max_signal_slope'] for r in rr],ls,label=f'{label}, out {ell}')
            axes[0,1].semilogy(range(1,51),[r['detuning_floor'] for r in rr],ls,label=f'{label}, out {ell}')
        axes[1,0].semilogy(range(1,51),state['fluctuation_occupancy'],label=label)
        for check in state['nonlinear']:
            rr=sorted([r for r in check['records'] if r['delta']>0],key=lambda r:r['delta'])
            if check['output_site']==50:
                axes[1,1].loglog([r['delta'] for r in rr],[r['relative_readout_error'] if r['converged'] else np.nan for r in rr],label=f"{label}, signal {check['signal_site']}, out 50")
    axes[0,0].set(xlabel='detuning signal site j',ylabel='maximum quadrature response / detuning')
    axes[0,1].set(xlabel='detuning signal site j',ylabel='optimal equivalent detuning ASD [sqrt(J)]')
    axes[1,0].set(xlabel='site',ylabel='linearized fluctuation photon number')
    axes[1,1].axhline(.05,color='k',ls=':');axes[1,1].set(xlabel='positive detuning shift / J',ylabel='relative nonlinear readout error')
    for ax in axes.flat:ax.legend(fontsize=7)
    fig.savefig('figures/compensated_sensing_pilot.png',dpi=180)


def report(s):
    states=s['states'];unc=states['uncompensated'];comp=states['compensated']
    lines=[r'''# Local-detuning sensing with and without boundary compensation

## Main finding

At N=50, epsilon=32, phi=pi/3, Delta=0.5, J=kappa=1, U=-0.0002 and V=0, the compensated intermediate state has a very large static signal response. **The amplified response is accompanied by amplified vacuum noise: large gain does not translate into a comparable improvement in detuning sensitivity.** More seriously, the linearized calculation predicts fluctuations much larger than the coherent background at the amplifying end. The compensated-state quantum sensitivity estimates below are therefore formal fixed-background BdG results, not a controlled prediction for a realizable quantum sensor at this operating point.

This study compares the previously prepared compensated homogeneous state (bulk winding -1 at omega=0) with the uncompensated vacuum attractor (trivial bulk). It uses mean-field backgrounds followed by Bogoliubov response/noise calculations. No fluctuation corrections are fed into the mean field and no covariance is evolved to find a stationary state.

## Measurement assumptions and definitions

Each site has an independent Markovian vacuum loss channel of rate kappa. Only one endpoint channel is measured at a time, with ideal, unit-efficiency homodyne detection; other output channels are unobserved. Coherent pumps are noiseless prescribed Hamiltonian drives. They are held fixed when a signal is applied, including the two boundary compensation fields. Any coherent drive background in a physical output is assumed calibrated/subtracted. We quote the coherent *emission contribution* kappa*|alpha|^2, not a total flux including interference with an incident pump. No detector noise, technical pump noise, fluctuating signal source or signal backaction is modeled. We estimate one known-position static detuning shift at a time, with an independently calibrated homodyne angle.

Define quadratures x=(delta a+delta a^dagger)/sqrt(2), p=(delta a-delta a^dagger)/(i*sqrt(2)), ordered r=(x_1,...,x_N,p_1,...,p_N). The full mean-field/BdG conventions remain i*dot(z)=L*z, with onsite U*n^2 and its +U shift retained. If A and B are the particle and pairing blocks, the real evolution generator is

```text
G = [[ Im(A+B),  Re(A-B)],
     [-Re(A+B),  Im(A-B)]].
```

For a local detuning shift deltaDelta_j, the signal force is f_j=sqrt(2)*(Im(alpha_j)*e_xj - Re(alpha_j)*e_pj). The zero-frequency intracavity response is R_j=-G^(-1)f_j. The Langevin equation and input-output relation are

```text
dot(r) = G*r + f_j*deltaDelta_j + sqrt(kappa)*r_in,
r_out  = r_in - sqrt(kappa)*r.
```

Thus the output signal slope is d_j=-sqrt(kappa)*R_j, while the zero-frequency vacuum scattering matrix is M=I+kappa*G^(-1). The full two-sided symmetrized output noise PSD matrix is C_out(0)=M*M^T/2. Vacuum quadrature PSD equals 1/2 in this convention. Keeping the identity in M includes the interference/correlations between directly reflected vacuum and cavity-emitted noise; simply adding shot noise to intracavity noise would be incorrect. These input-output and symmetrized-noise conventions follow the framework reviewed by [Clerk et al., Rev. Mod. Phys. 82, 1155 (2010)](https://arxiv.org/abs/0810.4729); the model-specific response matrices and numerical results here are new calculations.

At output site ell, let C_ell be its 2x2 output covariance spectral matrix and d_ell,j its two-component signal slope. A homodyne quadrature u=(cos(theta),sin(theta)) has slope u^T*d and PSD u^T*C*u. Its equivalent detuning amplitude spectral density is

```text
D_Delta = sqrt(u^T*C*u) / |u^T*d|,
u_opt proportional to C^(-1)*d,
D_Delta,opt = 1/sqrt(d^T*C^(-1)*d).
```

This noise-optimal quadrature generally differs from the maximum-signal quadrature u_gain=d/||d||. D_Delta has units sqrt(J). For a long averaging time T, the formal white-noise standard error is D_Delta/sqrt(T), with T in units 1/J. This is an ideal homodyne imprecision, not a quantum Fisher-information bound and not an optimization over joint multiport measurements.

## Signal and noise results

The signal site is scanned over all 50 sites, measuring each end separately and optimizing its homodyne angle. The best pair for each background is:

| Background | Signal site | Output site | theta (rad) | Signal slope | Noise PSD | Equivalent detuning ASD / sqrt(J) |
|:---|---:|---:|---:|---:|---:|---:|''']
    for label,state in states.items():
        b=state['best'];lines.append(f"| {label} | {b['signal_site']} | {b['output_site']} | {b['theta']:.6f} | {b['slope']:.6g} | {b['noise_psd']:.6g} | {b['detuning_floor']:.7g} |")
    lines += [r'''
The formal best compensated noise floor is only about 1.20 times better than the best uncompensated floor, despite the much larger unrestricted resolvent gain. These optima use different signal locations; they are not an equal-location comparison.

For the compensated state's most amplified route, signal at site 2 and output at site 50, the maximum-signal slope is 5.27078e5 and the noise PSD in that quadrature is 3.20788e8. Their ratio gives an equivalent detuning ASD of 0.0339808*sqrt(J), worse than the compensated state's low-gain optimal measurement at site 1. The noise-optimal and maximum-signal quadratures on this amplified route give essentially the same floor.

At the same signal site 2, using output site 1, the noise-optimal ASD is 0.0412172*sqrt(J) without compensation and 0.0198126*sqrt(J) with compensation: a formal factor 2.08. For a central signal at site 25, the best endpoint ASD changes from 13.0818*sqrt(J) without compensation to 4.86568*sqrt(J) with compensation (output 50 for the latter). This illustrates the strong dependence on both perturbation position and readout choice. It is not a universal topological sensing advantage.

![Signal response, noise-equivalent sensitivity, fluctuation occupancy and nonlinear readout check](../figures/compensated_sensing_pilot.png)

## Resources and a crucial consistency check

| Quantity | Uncompensated | Compensated |
|:---|---:|---:|''']
    for label,key in [('Total coherent photon population','population'),('Total coherent emission / J','total_coherent_emission')]:
        lines.append(f"| {label} | {unc[key]:.7g} | {comp[key]:.7g} |")
    lines += [f"| Left endpoint coherent emission / J | {unc['edge_flux'][0]:.7g} | {comp['edge_flux'][0]:.7g} |",
              f"| Right endpoint coherent emission / J | {unc['edge_flux'][1]:.7g} | {comp['edge_flux'][1]:.7g} |",
              f"| Maximum linearized fluctuation photons at a site | {max(unc['fluctuation_occupancy']):.7g} | {max(comp['fluctuation_occupancy']):.7g} |",
              f"| Maximum fluctuation/coherent population ratio | {unc['max_fluctuation_to_coherent']:.7g} | {comp['max_fluctuation_to_coherent']:.7g} |",r'''
The compensated state uses 4.85 times more coherent photons. A descriptive population-normalized metric D_Delta*sqrt(total coherent population), lower being better, is about 1.84 times worse for its individually optimized measurement. This is a resource comparison at the two operating points, not a demonstrated equal-resource optimum or an accounting of pump power and local-oscillator resources.

For consistency only, we solve the stationary *linearized* Lyapunov equation G*V+V*G^T=-(kappa/2)*I. Its fluctuation photon population is n_fl,j=(V_xjxj+V_pjpj-1)/2. This is a diagnostic after fixing the mean field, not a self-consistent Gaussian stationary solver.

The compensated state has roughly 1.43e7 predicted fluctuation photons at site 50, compared with a coherent density of 2988.6: a ratio about 4800. The total predicted fluctuation population is 3.37e7. **The small-fluctuation expansion fails this consistency check.** Negative BdG eigenvalue growth rates guarantee stability of the linear equations but do not guarantee small vacuum-driven fluctuations in a non-normal amplifier. Consequently none of the apparent compensated-state quantum sensitivity advantages here should be treated as established physical performance. The uncompensated state passes this check much better: its largest local fluctuation/coherent ratio is 4.18e-4.

## Nonlinear stationary response checks

For the optimized signal location of each state and for site 25, we apply positive and negative shifts from 1e-10*J to 0.1*J in half-decade steps. Each perturbed stationary solve starts from the unperturbed state plus its full linear prediction, with pumps unchanged. We check the nonlinear output shift against the predicted slope at a *fixed* noise-optimal homodyne angle calibrated at zero signal; output noise is not recomputed at finite signal. The following are the last contiguous sampled magnitudes, starting at the first numerically resolved weak-signal point (absolute predicted quadrature shift above 1e-8), for which the solver converges, the root is BdG stable, and the readout error is below 5%:

| Background | Signal site | Output site | Negative-shift magnitude / J | Positive-shift magnitude / J |
|:---|---:|---:|---:|---:|''']
    for label,state in states.items():
        for check in state['nonlinear']:
            limits=[]
            for sign in [-1,1]:
                rr=sorted([r for r in check['records'] if sign*r['delta']>0 and abs(r['prediction'])>1e-8],key=lambda r:abs(r['delta']));good=[]
                for r in rr:
                    if r['converged'] and r['growth']<0 and r['relative_readout_error']<.05:good.append(abs(r['delta']))
                    else:break
                limits.append(f'{max(good):.5g}' if good else 'none')
            lines.append(f"| {label} | {check['signal_site']} | {check['output_site']} | "+' | '.join(limits)+' |')
    lines += [r'''
These are sampled local-root checks, not certified dynamical switching thresholds. A solver failure does not establish loss of the branch, and a converged root need not be the attractor reached by suddenly applying a finite signal. No finite-signal time-evolution or rare switching simulation is claimed. The asymmetry between signal signs reflects nonlinear response and branch/basin structure. Mean-field agreement at weak signals validates the deterministic susceptibility, not the validity of the compensated quantum-noise expansion.

## Verification, reproducibility and interpretation

The real-quadrature susceptibility is independently checked against -L^(-1)*diag(alpha,-alpha*) in particle/hole coordinates. A passive U=0 reference has exactly vacuum output PSD=I/2, checking the direct-input interference and normalization. The all-port zero-frequency scattering matrix preserves quadrature commutators to maximum absolute error 4.9e-8 in the compensated case and below 1e-15 in the uncompensated case; cancellation involves large amplified coefficients. Relative Lyapunov residuals are below 1.2e-14. The reference stationary residuals are below 1e-9; all individual perturbed residuals, convergence flags and growth rates are retained, including failed solves. Densities and both original vacuum/ramp preparation histories come from the existing validated epsilon=32 dataset.

The calculation demonstrates a classical signal amplifier and shows why topology/large susceptibility alone do not establish quantum-enhanced sensing. A useful next step is to look for a shorter chain or a less strongly amplified operating point where fluctuation populations remain small, then repeat the same physical-signal and noise comparison. Technical compensation noise and finite efficiency would also need evaluation before drawing experimental conclusions.

Reproduce with `python src/compensated_sensing_pilot.py`. Full per-site sensitivities, output covariance spectra at zero frequency, optimal homodyne angles, nonlinear root checks and fluctuation diagnostics are saved in `data/generated/compensated_sensing/sensing.json`. Source mean fields and preparation data: `data/generated/compensation_n50/preparation_30_35/epsilon_32.000000.json`. Only this separate report is written by the script; the previous cumulative report is unchanged.
''']
    Path('docs/compensated_topological_sensing_pilot.md').write_text('\n'.join(lines),encoding='utf-8')


if __name__=='__main__':main()
