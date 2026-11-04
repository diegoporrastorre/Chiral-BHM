# Why compensation sensitivity grows exponentially with chain length

## Model, target, and meaning of sensitivity

Use J=1, U=-0.0002, V=0, Delta=0.5, kappa=1, phi=pi/3, epsilon=38.51, N=10..100. These restore the original V=0 baseline; the later Delta=0, kappa=4 studies are not used. The onsite Hamiltonian is U n_j^2, so its normal-ordering +U shift remains explicit. Directed hopping is T[j+1,j]=J exp(+i phi), and sites run from 1 to N.

The intermediate bulk root is n0=1990.520824756707, alpha0=-36.367537560479-25.844207020991i. The compensation additions are hL=J exp(+i phi) alpha0=4.197971040602-44.417314910955i and hR=J exp(-i phi) alpha0=-40.565508601082+18.573107889965i. Both have magnitude 44.615253275. They are independent of N. The total endpoint pumps are epsilon+hL and epsilon+hR.

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

For every integer N from 10 to 100, the target is stationary (maximum residual 3.24e-14) and its OBC BdG operator is stable (largest computed growth -0.051460 J). The fitted large-N response is

`S_common(N) ~= 1.305994 exp(0.167987171 N)` for N=50..100,

with R^2=0.999998802. The worst physical endpoint gain has exponent 0.167984878; inverse smallest singular value has exponent 0.167984880. The fitted slopes change little when the lower cutoff is moved from N=30 to 50, while N=10..30 shows visible oscillatory finite-size corrections.

|quantity|fit N interval|exponent per site|e-fold length (sites)|R squared|
|---|---|---:|---:|---:|
|common_relative_peak|10-100|0.168389819|5.938601|0.999742679|
|common_relative_peak|30-100|0.167986818|5.952848|0.999985159|
|common_relative_peak|50-100|0.167987171|5.952836|0.999998802|
|worst_endpoint_relative_peak|10-100|0.167619408|5.965896|0.999673236|
|worst_endpoint_relative_peak|30-100|0.167949470|5.954172|0.999984902|
|worst_endpoint_relative_peak|50-100|0.167984878|5.952917|0.999998802|
|inverse_smin|10-100|0.167721642|5.962260|0.999681211|
|inverse_smin|30-100|0.167949889|5.954157|0.999984906|
|inverse_smin|50-100|0.167984880|5.952917|0.999998802|

The common-error tolerance in the following table is the **linear estimate** eta_1%=0.01/S_common for keeping maximum complex-amplitude error at 1%. It is not a guaranteed basin size, a hard failure threshold, or an experimentally validated precision budget.

|N|S_common|linear common fractional error for 1% amplitude change|OBC growth / J|
|---:|---:|---:|---:|
|10|7.63402|0.00130993|-0.102504|
|20|33.771|0.000296112|-0.073533|
|30|203.661|4.91013e-05|-0.063741|
|40|1082.61|9.23695e-06|-0.059914|
|50|5780.41|1.72998e-06|-0.057127|
|60|31223.7|3.20269e-07|-0.055166|
|70|166804|5.99505e-08|-0.053971|
|80|896852|1.11501e-08|-0.052856|
|90|4.80663e+06|2.08046e-09|-0.052156|
|100|2.57976e+07|3.87633e-10|-0.051460|

At N=50 the estimate is 1.72998e-06, about 1.7300 ppm; at N=100 it is 3.87633e-10, about 0.3876 ppb. The latter corresponds to an absolute error in either compensation addition of approximately 1.72943e-08 J. The required nominal compensation itself has not grown with N.

### Direction and quadrature matter

|N|left amplitude gain|left phase gain (per radian)|right amplitude gain|right phase gain (per radian)|worst unit-norm endpoint gain|
|---:|---:|---:|---:|---:|---:|
|10|10.1495|12.3028|2.51728|0.499482|16.1463|
|50|5781.58|7008.8|1.1708|0.540305|9085.7|
|100|2.57976e+07|3.12735e+07|1.16467|0.546346|4.05407e+07|

For this hopping orientation, left-end errors are amplified to the right. Right-end gains approach constants of order one. The common two-end amplitude error is asymptotically dominated by its left contribution. Reversing the chirality and reflecting the chain would exchange the roles. Thus it would be misleading to assign the same exponential sensitivity to every endpoint and every input direction. A special input quadrature with zero overlap with the amplified channel need not share its leading exponential factor.

## Spatial-root explanation

After a global Nambu rotation the constant pairing can be chosen real, |g|=|2U|n0. Let d=Delta+U+4Un0 and gamma=kappa/2. The bulk symbol for Fourier convention exp(ikj) is

`L(k)=(2J sin(phi) sin(k)-i gamma)I+[[d+2J cos(phi)cos(k),g],[-g,-d-2J cos(phi)cos(k)]]`.

A zero-frequency stationary bulk mode has dependence beta^j, beta=exp(ik), and det L(beta)=0. Equivalently, solve the lossless symbol at omega=i gamma. Multiplying the determinant by beta^2 gives a quartic, with roots

|beta|abs(beta)|
|---|---:|
|0.504034057 +0.000000000i|0.504034057|
|1.182920770 +0.000000000i|1.182920770|
|-0.730181786 -1.069594803i|1.295066980|
|-0.730181786 +1.069594803i|1.295066980|

There is one root inside and three outside the unit circle, consistent with zero-frequency winding -1 in this orientation. Let beta_star be the outside root closest to the unit circle:

`beta_star=1.182920770275`, `mu=log|beta_star|=0.167986609191`, `xi=1/mu=5.952855438 sites`.

Here is the boundary-mode argument for that length scale. Combine the three outside-root modes as psi_j=sum_(l=2..4) c_l chi_l beta_l^(j-(N+1)), choosing their coefficients to satisfy psi_(N+1)=0. There are three coefficients and two right-boundary conditions, leaving a nonzero right-localized vector. Its bulk recurrence and right boundary are satisfied exactly. After truncating to sites 1..N, the residual at the left edge is controlled generically by the slowest tail, |beta_star|^(-N). Therefore it constructs a quasi-zero singular vector with s_min bounded above by a constant times exp(-mu N). Provided the leading boundary overlap does not vanish, this gives the generic asymptotic s_min~exp(-mu N), and physical endpoint errors overlapping its left singular vector have gain~exp(+mu N).

For this parameter set the measured singular gap and physical endpoint responses confirm that asymptotic scaling, including the same exponent. The argument identifies the controlling spatial root; it is not a theorem that every observable, input direction, or arbitrary inhomogeneous chain has this exponent. The other two outside roots have larger modulus and produce faster-decaying, oscillatory finite-size corrections. Their relative magnitude predicts a correction envelope approximately (|beta_star|/1.29506698)^N, explaining why small-N oscillations diminish.

Each added site asymptotically multiplies the amplified response by 1.182921; every 13.706956 sites adds roughly one decimal order of sensitivity. The e-fold length is 5.952855 sites. These are **spatial** amplification/attenuation scales, not temporal instability growth rates.

## Why this is compatible with stability

All tested OBC growth rates remain negative. Local dynamical stability concerns poles/eigenvalues, whereas static sensitivity concerns the resolvent at omega=0 and singular values. A nonnormal operator can have decaying eigenmodes and an exponentially small singular value. Consequently the precisely compensated state is locally attracting, yet a persistent boundary-drive mismatch can displace its stationary background enormously. Linear stability alone does not guarantee robustness, a large basin, or successful vacuum preparation.

## Independent precision and finite-error verification

At N=10,20,40,50,60,80,100 solve the physical four-channel stationary derivative independently with a real banded Decimal Jacobian at 100 and 160 digits, recomputing the bulk root and compensation at high precision. Both precisions agree; the largest relative discrepancy between double precision and the 160-digit common gain is 1.15e-13. These checks verify the static response directly rather than relying only on an ill-conditioned double-precision inverse or SVD.

At N=10,50,100, independently solve the **nonlinear stationary equations** with both signs of a common compensation error, choosing eta to predict relative amplitude changes 1e-5 and 1e-3. High-precision banded Newton iterations reach force residuals below 1e-55. The largest fractional deviation of actual response from the linear prediction is 0.00156647, about 0.1566 percent. This supports the derivative in its small-error regime, including the extremely small applied eta at N=100. It does not establish the 1% error-budget estimate as an exact nonlinear threshold, and these checks are stationary solves, not new time-evolution or noise simulations.

Code: src/compensation_sensitivity_scaling.py and src/compensation_precision.py. Dataset: data/generated/compensation_sensitivity/scaling.json. Run `python src/compensation_sensitivity_scaling.py`. Figures are figures/compensation_sensitivity_scaling.png and figures/compensation_sensitivity_exponent.png. Parameter signs, input normalizations, fit windows, residuals, growth rates, spatial roots, precision comparisons and nonlinear finite-error checks are recorded with the results. No covariance evolution or fluctuation feedback is included.
