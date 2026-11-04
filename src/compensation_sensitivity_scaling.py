"""N=10..100: static sensitivity to physical endpoint compensation errors."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json
import numpy as np
from scipy.linalg import solve,svdvals,eigvals
from mean_field import Chain
from bulk_stability import spatial_roots
from compensation_precision import derivatives,nonlinear_check
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DEST=Path('data/generated/compensation_sensitivity');DEST.mkdir(parents=True,exist_ok=True)

def fit(rows,key,minimum):
    rr=[r for r in rows if r['N']>=minimum];x=np.array([r['N'] for r in rr]);y=np.log([r[key] for r in rr])
    slope,intercept=np.polyfit(x,y,1);res=y-(intercept+slope*x)
    return dict(key=key,min_N=minimum,max_N=int(x.max()),slope=float(slope),prefactor=float(np.exp(intercept)),
                length=float(1/slope),r_squared=float(1-sum(res*res)/sum((y-y.mean())**2)),
                maximum_log_residual=float(max(abs(res))))

def main():
    epsilon=38.51;c=1.4998;beta=-.0004
    n=sorted(np.roots([beta*beta,2*c*beta,c*c+.25,-epsilon**2]).real)[1]
    alpha=-epsilon/(c+beta*n-.5j);phi=np.pi/3;hl=np.exp(1j*phi)*alpha;hr=np.exp(-1j*phi)*alpha
    rows=[]
    for N in range(10,101):
        model=Chain(N=N,drive=epsilon);model.pump=model.pump.astype(complex);model.pump[0]+=hl;model.pump[-1]+=hr
        a=np.full(N,alpha);residual=float(max(abs(model.force(a))));assert residual<1e-9
        A,B=model.blocks(a)
        K=np.block([[(A+B).real,-(A-B).imag],[(A+B).imag,(A-B).real]])
        q=np.zeros((N,4),complex);q[0,0]=hl;q[0,1]=1j*hl;q[-1,2]=hr;q[-1,3]=1j*hr
        X=solve(K,-np.r_[q.real,q.imag]);R=X[:N]+1j*X[N:];common=R[:,0]+R[:,2]
        channel=np.max(abs(R),axis=0)/abs(alpha)
        worst=max(np.linalg.svd(np.vstack([R[j].real,R[j].imag]),compute_uv=False)[0] for j in range(N))/abs(alpha)
        density_common=2*np.real(alpha.conjugate()*common)/n
        L=model.bogoliubov(a);singular=svdvals(L)
        row=dict(N=N,residual=residual,growth=float(eigvals(L).imag.max()),smin=float(singular[-1]),inverse_smin=float(1/singular[-1]),
                 common_peak=float(max(abs(common))),common_relative_peak=float(max(abs(common)))/abs(alpha),
                 common_relative_density_peak=float(max(abs(density_common))),
                 left_amplitude=float(channel[0]),left_phase=float(channel[1]),right_amplitude=float(channel[2]),right_phase=float(channel[3]),
                 worst_endpoint_relative_peak=float(worst),
                 solve_relative_residual=float(np.linalg.norm(K@X+np.r_[q.real,q.imag])/np.linalg.norm(np.r_[q.real,q.imag])),
                 common_relative_profile=(abs(common)/abs(alpha)).tolist())
        row['common_error_for_one_percent_amplitude']=.01/row['common_relative_peak']
        row['worst_error_for_one_percent_amplitude']=.01/worst
        rows.append(row)
    d=.4998-.0008*n;g=.0004*n
    roots=sorted(spatial_roots(.5j,1,phi,d,g),key=abs)
    selected_root=min([z for z in roots if abs(z)>1],key=abs);mu=float(np.log(abs(selected_root)))
    result=dict(parameters=dict(J=1,U=-.0002,V=0,delta=.5,kappa=1,epsilon=epsilon,phi=float(phi)),
                density=float(n),alpha0=[float(alpha.real),float(alpha.imag)],
                compensation_left=[float(hl.real),float(hl.imag)],compensation_right=[float(hr.real),float(hr.imag)],
                rows=rows,spatial_roots=[[float(z.real),float(z.imag)] for z in roots],predicted_exponent=mu,predicted_length=1/mu,
                fits=[fit(rows,key,minimum) for key in ['common_relative_peak','worst_endpoint_relative_peak','inverse_smin'] for minimum in [10,30,50]],
                precision_checks=[],nonlinear_checks=[])
    for N in [10,20,40,50,60,80,100]:
        high=[derivatives(N,p) for p in [100,160]]
        assert abs(high[0]['common_relative_peak']/high[1]['common_relative_peak']-1)<1e-13
        row=next(r for r in rows if r['N']==N)
        diff=abs(row['common_relative_peak']/high[1]['common_relative_peak']-1)
        result['precision_checks'].append(dict(N=N,double_relative_difference=float(diff),checks=high))
        print('precision',N,high[-1]['common_relative_peak'],diff,flush=True)
    for N in [10,50,100]:
        sensitivity=next(r for r in result['precision_checks'] if r['N']==N)['checks'][-1]['common_relative_peak']
        for target in [1e-5,1e-3]:
            for sign in [-1,1]:
                eta=sign*target/sensitivity;check=nonlinear_check(N,eta,100)
                check.update(linear_prediction=target,actual_to_linear_ratio=check['relative_amplitude_change']/target)
                result['nonlinear_checks'].append(check)
    (DEST/'scaling.json').write_text(json.dumps(result,indent=2))
    plot(result);report(result)
    print('fits',result['fits'],'predicted exponent',mu,flush=True)

def plot(result):
    rr=result['rows'];N=np.array([r['N'] for r in rr]);mu=result['predicted_exponent']
    fig,ax=plt.subplots(2,2,figsize=(12,9),layout='constrained')
    for key,label in [('common_relative_peak','common amplitude error'),('left_amplitude','left amplitude error'),('left_phase','left phase error'),('right_amplitude','right amplitude error'),('right_phase','right phase error')]:
        ax[0,0].semilogy(N,[r[key] for r in rr],label=label)
    f=next(f for f in result['fits'] if f['key']=='common_relative_peak' and f['min_N']==50)
    ax[0,0].semilogy(N,f['prefactor']*np.exp(f['slope']*N),'k--',label='common error: N>=50 exponential fit')
    ax[0,0].set(xlabel='N',ylabel='peak relative amplitude response / error');ax[0,0].legend(fontsize=8)
    ax[0,1].semilogy(N,[r['inverse_smin'] for r in rr],label='1 / smallest singular value (J units)')
    ax[0,1].semilogy(N,[r['worst_endpoint_relative_peak'] for r in rr],label='worst physical endpoint-error gain');ax[0,1].legend(fontsize=8);ax[0,1].set(xlabel='N',ylabel='static gain')
    ax[1,0].semilogy(N,[r['common_error_for_one_percent_amplitude'] for r in rr],label='common fractional amplitude error')
    ax[1,0].semilogy(N,[r['worst_error_for_one_percent_amplitude'] for r in rr],label='worst unit-norm endpoint error');ax[1,0].legend(fontsize=8);ax[1,0].set(xlabel='N',ylabel='linear error budget for 1% amplitude change')
    for size in [10,20,50,100]:
        r=next(r for r in rr if r['N']==size);ax[1,1].semilogy(range(1,size+1),r['common_relative_profile'],label=f'N={size}')
    ax[1,1].set(xlabel='site j',ylabel='relative response to common error');ax[1,1].legend()
    for axes in ax.flat:axes.grid(alpha=.15)
    fig.savefig('figures/compensation_sensitivity_scaling.png',dpi=180)
    fig,ax=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    for key,label in [('common_relative_peak','common-error gain'),('worst_endpoint_relative_peak','worst endpoint gain'),('inverse_smin','inverse singular gap')]:
        y=np.array([r[key] for r in rr]);ax[0].plot(N,y*np.exp(-mu*N),label=label);ax[1].plot((N[1:]+N[:-1])/2,np.diff(np.log(y)),label=label)
    ax[0].set(xlabel='N',ylabel='gain divided by exp(mu N)');ax[1].axhline(mu,color='k',ls='--',label='spatial-root prediction');ax[1].set(xlabel='N',ylabel='successive log-gain increment');ax[0].legend(fontsize=8);ax[1].legend(fontsize=8)
    fig.savefig('figures/compensation_sensitivity_exponent.png',dpi=180)
    fig,ax=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    fig.suptitle(r'Compensation sensitivity: $V=0$, $\Delta=0.5$, $\kappa=J=1$, $U=-2\times10^{-4}$'+'\n'+r'$\phi=\pi/3$, $\epsilon=38.51$; homogeneous intermediate state')
    ax[0].semilogy(N,[r['common_relative_peak'] for r in rr],'-o',ms=3,markevery=5,label='stationary derivative')
    ax[0].semilogy(N,f['prefactor']*np.exp(f['slope']*N),'k--',label=rf"fit: $1.306e^{{{f['slope']:.5f}N}}$")
    ax[0].set(xlabel='Chain length N',ylabel=r'Common-error sensitivity $S_N$');ax[0].legend(fontsize=9)
    ax[1].semilogy(N,[r['common_error_for_one_percent_amplitude'] for r in rr],'-o',ms=3,markevery=5)
    ax[1].set(xlabel='Chain length N',ylabel='Fractional compensation-error budget',title='Linear estimate for 1% amplitude change')
    for axes in ax:axes.grid(alpha=.2)
    fig.savefig('figures/compensation_sensitivity_overview.png',dpi=180)

def report(result):
    rows=result['rows'];common_fit=next(f for f in result['fits'] if f['key']=='common_relative_peak' and f['min_N']==50)
    worst_fit=next(f for f in result['fits'] if f['key']=='worst_endpoint_relative_peak' and f['min_N']==50)
    singular_fit=next(f for f in result['fits'] if f['key']=='inverse_smin' and f['min_N']==50)
    mu=result['predicted_exponent'];xi=result['predicted_length'];alpha=complex(*result['alpha0'])
    hl=complex(*result['compensation_left']);hr=complex(*result['compensation_right'])
    sample=[]
    for r in rows:
        if r['N']%10==0:
            sample.append(f"|{r['N']}|{r['common_relative_peak']:.6g}|{r['common_error_for_one_percent_amplitude']:.6g}|{r['growth']:.6f}|")
    channel=[]
    for N in [10,50,100]:
        r=next(r for r in rows if r['N']==N)
        channel.append(f"|{N}|{r['left_amplitude']:.6g}|{r['left_phase']:.6g}|{r['right_amplitude']:.6g}|{r['right_phase']:.6g}|{r['worst_endpoint_relative_peak']:.6g}|")
    fit_rows=[f"|{f['key']}|{f['min_N']}-{f['max_N']}|{f['slope']:.9f}|{f['length']:.6f}|{f['r_squared']:.9f}|" for f in result['fits']]
    r100=rows[-1];r50=next(r for r in rows if r['N']==50)
    high_difference=max(z['double_relative_difference'] for z in result['precision_checks'])
    nonlinear_error=max(abs(z['actual_to_linear_ratio']-1) for z in result['nonlinear_checks'])
    roots=[complex(*z) for z in result['spatial_roots']]
    roots_table=[f"|{z.real:.9f} {z.imag:+.9f}i|{abs(z):.9f}|" for z in roots]
    derivation=f'''# Why compensation sensitivity grows exponentially with chain length

## Model, target, and meaning of sensitivity

Use J=1, U=-0.0002, V=0, Delta=0.5, kappa=1, phi=pi/3, epsilon=38.51, N=10..100. These restore the original V=0 baseline; the later Delta=0, kappa=4 studies are not used. The onsite Hamiltonian is U n_j^2, so its normal-ordering +U shift remains explicit. Directed hopping is T[j+1,j]=J exp(+i phi), and sites run from 1 to N.

The intermediate bulk root is n0={result['density']:.12f}, alpha0={alpha.real:.12f}{alpha.imag:+.12f}i. The compensation additions are hL=J exp(+i phi) alpha0={hl.real:.12f}{hl.imag:+.12f}i and hR=J exp(-i phi) alpha0={hr.real:.12f}{hr.imag:+.12f}i. Both have magnitude {abs(alpha):.9f}. They are independent of N. The total endpoint pumps are epsilon+hL and epsilon+hR.

The N-dependent question is how much the stationary profile changes if those fixed additions are imperfect. Let actual additions be (1+eta)hL and (1+eta)hR for a common real fractional amplitude error eta. Define

`S_common(N)=max_j |d alpha_j/d eta| / |alpha0|`.

Then max_j|delta alpha_j|/|alpha0|=S_common|eta|+O(eta^2). This measures complex-amplitude distance, including phase changes; it is not just a density error or a tolerance on the total endpoint pump. For density, d n_j/d eta=2 Re(alpha0* d alpha_j/d eta).

## Exact stationary derivative

For the Nambu perturbation z=(delta alpha,delta alpha*) and forcing q=d pump/d eta, differentiating the stationary equations gives

`L z' = -(q,-q*)`,

where q_1=hL, q_N=hR, other q_j=0, and

`L=[[A,B],[-B*,-A*]]`,

`A=T+(Delta+U+4U n0-i kappa/2) I`, `B=2U alpha0^2 I`.

The equivalent real stationary Jacobian is K=[[Re(A+B),-Im(A-B)],[Im(A+B),Re(A-B)]]. Solve K (Re alpha',Im alpha')^T=-(Re q,Im q)^T. This fixes the conjugation and forcing signs and uses physical coherent drive errors; particle/hole inputs are not independently varied.

For independent small errors, write delta hL=hL(etaL+i thetaL), delta hR=hR(etaR+i thetaR). Theta is in radians and eta is fractional field amplitude. A real 2N-by-4 response matrix R maps x=(etaL,thetaL,etaR,thetaR) to the stationary deformation. The worst endpoint gain is max_j sigma_max(R_j)/|alpha0|, where R_j is the real 2-by-4 block at site j. Thus its input normalization is ||x||_2=1. The common-amplitude protocol uses x=(eta,0,eta,0), whose norm is sqrt(2)|eta|; keep these normalizations distinct.

## Numerical scaling

For every integer N from 10 to 100, the target is stationary (maximum residual {max(r['residual'] for r in rows):.3g}) and its OBC BdG operator is stable (largest computed growth {max(r['growth'] for r in rows):.6f} J). The fitted large-N response is

`S_common(N) ~= {common_fit['prefactor']:.6f} exp({common_fit['slope']:.9f} N)` for N=50..100,

with R^2={common_fit['r_squared']:.9f}. The worst physical endpoint gain has exponent {worst_fit['slope']:.9f}; inverse smallest singular value has exponent {singular_fit['slope']:.9f}. The fitted slopes change little when the lower cutoff is moved from N=30 to 50, while N=10..30 shows visible oscillatory finite-size corrections.

|quantity|fit N interval|exponent per site|e-fold length (sites)|R squared|
|---|---|---:|---:|---:|
{chr(10).join(fit_rows)}

The common-error tolerance in the following table is the **linear estimate** eta_1%=0.01/S_common for keeping maximum complex-amplitude error at 1%. It is not a guaranteed basin size, a hard failure threshold, or an experimentally validated precision budget.

|N|S_common|linear common fractional error for 1% amplitude change|OBC growth / J|
|---:|---:|---:|---:|
{chr(10).join(sample)}

At N=50 the estimate is {r50['common_error_for_one_percent_amplitude']:.6g}, about {r50['common_error_for_one_percent_amplitude']*1e6:.4f} ppm; at N=100 it is {r100['common_error_for_one_percent_amplitude']:.6g}, about {r100['common_error_for_one_percent_amplitude']*1e9:.4f} ppb. The latter corresponds to an absolute error in either compensation addition of approximately {abs(alpha)*r100['common_error_for_one_percent_amplitude']:.6g} J. The required nominal compensation itself has not grown with N.

### Direction and quadrature matter

|N|left amplitude gain|left phase gain (per radian)|right amplitude gain|right phase gain (per radian)|worst unit-norm endpoint gain|
|---:|---:|---:|---:|---:|---:|
{chr(10).join(channel)}

For this hopping orientation, left-end errors are amplified to the right. Right-end gains approach constants of order one. The common two-end amplitude error is asymptotically dominated by its left contribution. Reversing the chirality and reflecting the chain would exchange the roles. Thus it would be misleading to assign the same exponential sensitivity to every endpoint and every input direction. A special input quadrature with zero overlap with the amplified channel need not share its leading exponential factor.

## Spatial-root explanation

After a global Nambu rotation the constant pairing can be chosen real, |g|=|2U|n0. Let d=Delta+U+4Un0 and gamma=kappa/2. The bulk symbol for Fourier convention exp(ikj) is

`L(k)=(2J sin(phi) sin(k)-i gamma)I+[[d+2J cos(phi)cos(k),g],[-g,-d-2J cos(phi)cos(k)]]`.

A zero-frequency stationary bulk mode has dependence beta^j, beta=exp(ik), and det L(beta)=0. Equivalently, solve the lossless symbol at omega=i gamma. Multiplying the determinant by beta^2 gives a quartic, with roots

|beta|abs(beta)|
|---|---:|
{chr(10).join(roots_table)}

There is one root inside and three outside the unit circle, consistent with zero-frequency winding -1 in this orientation. Let beta_star be the outside root closest to the unit circle:

`beta_star={abs(roots[1]):.12f}`, `mu=log|beta_star|={mu:.12f}`, `xi=1/mu={xi:.9f} sites`.

Here is the boundary-mode argument for that length scale. Combine the three outside-root modes as psi_j=sum_(l=2..4) c_l chi_l beta_l^(j-(N+1)), choosing their coefficients to satisfy psi_(N+1)=0. There are three coefficients and two right-boundary conditions, leaving a nonzero right-localized vector. Its bulk recurrence and right boundary are satisfied exactly. After truncating to sites 1..N, the residual at the left edge is controlled generically by the slowest tail, |beta_star|^(-N). Therefore it constructs a quasi-zero singular vector with s_min bounded above by a constant times exp(-mu N). Provided the leading boundary overlap does not vanish, this gives the generic asymptotic s_min~exp(-mu N), and physical endpoint errors overlapping its left singular vector have gain~exp(+mu N).

For this parameter set the measured singular gap and physical endpoint responses confirm that asymptotic scaling, including the same exponent. The argument identifies the controlling spatial root; it is not a theorem that every observable, input direction, or arbitrary inhomogeneous chain has this exponent. The other two outside roots have larger modulus and produce faster-decaying, oscillatory finite-size corrections. Their relative magnitude predicts a correction envelope approximately (|beta_star|/1.29506698)^N, explaining why small-N oscillations diminish.

Each added site asymptotically multiplies the amplified response by {abs(roots[1]):.6f}; every {np.log(10)/mu:.6f} sites adds roughly one decimal order of sensitivity. The e-fold length is {xi:.6f} sites. These are **spatial** amplification/attenuation scales, not temporal instability growth rates.

## Why this is compatible with stability

All tested OBC growth rates remain negative. Local dynamical stability concerns poles/eigenvalues, whereas static sensitivity concerns the resolvent at omega=0 and singular values. A nonnormal operator can have decaying eigenmodes and an exponentially small singular value. Consequently the precisely compensated state is locally attracting, yet a persistent boundary-drive mismatch can displace its stationary background enormously. Linear stability alone does not guarantee robustness, a large basin, or successful vacuum preparation.

## Independent precision and finite-error verification

At N=10,20,40,50,60,80,100 solve the physical four-channel stationary derivative independently with a real banded Decimal Jacobian at 100 and 160 digits, recomputing the bulk root and compensation at high precision. Both precisions agree; the largest relative discrepancy between double precision and the 160-digit common gain is {high_difference:.3g}. These checks verify the static response directly rather than relying only on an ill-conditioned double-precision inverse or SVD.

At N=10,50,100, independently solve the **nonlinear stationary equations** with both signs of a common compensation error, choosing eta to predict relative amplitude changes 1e-5 and 1e-3. High-precision banded Newton iterations reach force residuals below 1e-55. The largest fractional deviation of actual response from the linear prediction is {nonlinear_error:.6g}, about {100*nonlinear_error:.4f} percent. This supports the derivative in its small-error regime, including the extremely small applied eta at N=100. It does not establish the 1% error-budget estimate as an exact nonlinear threshold, and these checks are stationary solves, not new time-evolution or noise simulations.

Code: src/compensation_sensitivity_scaling.py and src/compensation_precision.py. Dataset: data/generated/compensation_sensitivity/scaling.json. Run `python src/compensation_sensitivity_scaling.py`. Figures are figures/compensation_sensitivity_scaling.png and figures/compensation_sensitivity_exponent.png. Parameter signs, input normalizations, fit windows, residuals, growth rates, spatial roots, precision comparisons and nonlinear finite-error checks are recorded with the results. No covariance evolution or fluctuation feedback is included.
'''
    Path('derivations/compensation_sensitivity_scaling.md').write_text(derivation,encoding='utf-8')
    summary=f'''## Compensation sensitivity versus chain length: V=0, N=10..100

Return to **Delta=0.5, kappa=J=1, U=-0.0002, V=0**, with phi=pi/3 and epsilon=38.51. These are the original onsite-only baseline parameters, not the subsequent Delta=0/kappa=4 case. The intermediate density is n0={result['density']:.9f}; the nominal compensations hL=J exp(+i phi) alpha0 and hR=J exp(-i phi) alpha0 have magnitude {abs(alpha):.6f} and are **independent of N**.

What scales exponentially is the stationary response to an error in those additions. For common fractional amplitude error eta, pump_1=epsilon+(1+eta)hL and pump_N=epsilon+(1+eta)hR, define S_common=max_j|d alpha_j/d eta|/|alpha0|. Calculating every integer N=10..100 gives

**S_common(N) approximately {common_fit['prefactor']:.6f} exp({common_fit['slope']:.9f} N)** for N=50..100, with R^2={common_fit['r_squared']:.9f}.

![Exponential compensation sensitivity and the corresponding linear precision requirement versus chain length](../figures/compensation_sensitivity_overview.png)

![Compensation-error sensitivity, directional amplitude and phase channels, singular response and linear error budgets versus N](../figures/compensation_sensitivity_scaling.png)

### Meaning of the error tolerance

To keep the peak complex-amplitude change at 1%, the linear estimate is |eta| approximately 0.01/S_common. This is a fractional error of the **compensation addition**, not of the total endpoint pump or density. It is a small-signal estimate, not a hard nonlinear failure threshold.

|N|S_common|linear common fractional error for 1% amplitude change|OBC growth / J|
|---:|---:|---:|---:|
{chr(10).join(sample)}

The N=50 budget is about {r50['common_error_for_one_percent_amplitude']*1e6:.4f} ppm; N=100 is about {r100['common_error_for_one_percent_amplitude']*1e9:.4f} ppb. For the chosen orientation the growing sensitivity is primarily to **left-end** amplitude and phase errors. The right-end error responses remain order one. Independent phase-error parameters are in radians; the worst gain uses unit Euclidean norm of (etaL,thetaL,etaR,thetaR), and therefore has a different normalization from the common two-end amplitude protocol.

### Origin of the exponent

Differentiating the stationary equations gives L z'=-(q,-q*), so the response is a physical endpoint projection of the zero-frequency inverse BdG operator. Its smallest singular value decreases exponentially even though every tested OBC eigenmode decays (largest growth {max(r['growth'] for r in rows):.6f} J).

The zero-frequency spatial roots beta of the bulk BdG recurrence are approximately 0.50403406, 1.18292077 and -0.73018179 +/- 1.06959480i. The three outside-unit-circle modes combine into a right-localized quasi-zero singular vector. Its longest tail is set by beta_star=1.18292077, giving

**mu=log|beta_star|={mu:.9f} per site, xi=1/mu={xi:.6f} sites**.

Its boundary residual scales as exp(-mu N); endpoint inputs overlapping the corresponding left singular vector therefore have static gain exp(+mu N). The fitted physical response exponent {common_fit['slope']:.9f} agrees with that independent spatial-root prediction. One decimal order of extra sensitivity is added approximately every {np.log(10)/mu:.4f} sites. Oscillatory small-N corrections come from the faster-decaying complex spatial roots.

![Exponential scaling after removal of the predicted spatial factor, and convergence of successive log-gain increments](../figures/compensation_sensitivity_exponent.png)

### Verification and interpretation

Independent 100/160-digit banded derivative solves at N=10,20,40,50,60,80,100 agree with double precision to maximum relative difference {high_difference:.3g} for the common gain. Nonlinear high-precision stationary checks at N=10,50,100, both signs of error and predicted relative amplitude changes 1e-5 and 1e-3, agree with linear response within {100*nonlinear_error:.4f} percent. Force residuals are below 1e-55 in those nonlinear checks. The unperturbed compensated roots remain stable for all tested lengths.

The nominal compensation does not require an exponentially increasing field. Instead, the permitted **mismatch** shrinks exponentially if a fixed homogeneity tolerance is imposed on this amplified input channel. This is compatible with local dynamical stability and does not quantify global basins, time-dependent noise, or preparation success. Full input definitions, fit-window comparisons, spatial-root derivation, directional gains, and numerical checks are in [the analytical derivation](../derivations/compensation_sensitivity_scaling.md).

Reproduce with `python src/compensation_sensitivity_scaling.py`; helper src/compensation_precision.py implements the independent Decimal solves. Full records are data/generated/compensation_sensitivity/scaling.json. No covariance dynamics or fluctuation feedback is used.
'''
    path=Path('docs/uniform_topological_segment_stability.md');base=path.read_text(encoding='utf-8').split('## Compensation sensitivity versus chain length: V=0, N=10..100')[0].rstrip();path.write_text(base+'\n\n'+summary,encoding='utf-8')

if __name__=='__main__':main()
