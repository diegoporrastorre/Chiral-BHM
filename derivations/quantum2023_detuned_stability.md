# Coupled channels: Delta=0.1 in the Quantum 2023 convention

## Model and scope

This repeats `quantum2023_kappa_stability.md` with Delta changed from 0 to
0.1. All other conventions remain fixed: gs=gc=J=1, phi=pi/2 exactly,
collective pump P=0, local Lindblad loss kappa, amplitude loss kappa/2,
hbar=1, lattice spacing one. We use the dynamical-matrix convention of
Gómez-León et al., [Quantum 7, 1016 (2023)](https://doi.org/10.22331/q-2023-05-23-1016),
Eqs. (7)-(10) and (28)-(31), including the onsite-pairing normalization
discussed in the preceding derivation. No mean-field solution or Gaussian
covariance evolution is involved: this compares boundary conditions for
the same prescribed quadratic operator.

With u_j proportional to exp(ikj),

$$
H(k)=(-2\sin k-i\kappa/2)I+\delta\sigma_z+i(1+2\cos k)\sigma_y,
\qquad \delta=0.1.
$$

Evolution is i partial_t Psi=H Psi. Growth is max Im omega, equivalently
max Re lambda for A=-iH. Stability throughout means strict exponential
decay, not merely a nonpositive exponent. Sites are j=0,...,N-1. OBC
simply removes the closing bond; PBC uses the same phase on that bond,
with directed total phase -N phi in the paper's lower-diagonal convention.

The same constant unitary rotation that diagonalizes sigma_y now gives,
with z=exp(ik) and temporarily kappa=0,

$$
A(z)=\begin{pmatrix}1+2z&-i\delta\\-i\delta&-1-2/z\end{pmatrix}.
$$

Thus detuning converts one directed channel into the other. The finite
OBC matrix has z replaced by S, S[j,j+1]=1, and z^{-1} by S^T. Uniform
loss subtracts kappa/2 from every growth eigenvalue, so a loss sweep is
exactly a shift of one lossless calculation, not repeated ill-conditioned
diagonalizations at every loss.

## 1. PBC threshold

The Bloch frequencies are

$$
\omega_\pm(k)=-2\sin k-i\kappa/2
\pm\sqrt{\delta^2-(1+2\cos k)^2}.
$$

The largest imaginary part occurs at k=0 (present for any periodic N):

$$
\Gamma_{\rm PBC}=\sqrt{9-\delta^2}-\kappa/2,
\qquad \kappa_{c,\rm PBC}=2\sqrt{8.99}=5.9966657402.
$$

## 2. Local bulk response: finite saddles alone are incomplete

Write the characteristic Laurent polynomial for the lossless generator:

$$
D(\lambda,z)=\det[\lambda I-A(z)]
=a+bz+c/z,
$$

$$
a=\lambda^2-5+\delta^2,\quad b=-2(\lambda+1),\quad c=2(\lambda-1).
$$

Its spatial-root discriminant is

$$
F(\lambda)=a^2-4bc
=\lambda^4+(6+2\delta^2)\lambda^2+9-10\delta^2+\delta^4.
$$

The finite double roots obey lambda^2=-3-delta^2 +/-4delta, giving
lambda=+/-i sqrt(2.61), +/-i sqrt(3.41). Their real parts vanish.
However, it would be incorrect to infer zero local lossless growth from
these saddles alone: the spatial hopping is rank deficient, and additional
singularities appear when a spatial root reaches zero or infinity.

For Re lambda sufficiently large the on-site Laplace resolvent is

$$
R_{00}(\lambda)=\oint_{|z|=1}\frac{dz}{2\pi iz}
[\lambda I-A(z)]^{-1}.
$$

Choose s=sqrt(F) with s asymptotic to lambda^2+3+delta^2 at positive
large lambda, and continue from there. The contour identities are

$$
I_0=\oint\frac{dz}{2\pi izD}=\frac1s,\qquad
I_{-1}=\oint\frac{dz}{2\pi iz}\frac{z^{-1}}D
=\frac{1-a/s}{2c}.
$$

They follow from a I0+b I1+c I(-1)=1 and b I1=c I(-1), the latter
being the integral of a logarithmic derivative. Therefore

$$
R_{00,++}(\lambda)=\frac{\lambda^2+3-\delta^2+s}
{2(\lambda-1)s},\qquad
R_{00,+-}(\lambda)=-\frac{i\delta}{s},
$$

and R00,-- has the same numerator as R00,++ with lambda+1 in the
denominator. In the half-plane Re lambda>0 there is a pole at lambda=1
with residue one in R00,++; at delta=0.1 the branch points are all on
the imaginary axis. Thus the rightmost local singularity is lambda=1:

$$
\boxed{\Gamma_{\rm bulk}(v=0)=1-\kappa/2.}
$$

The local bulk criterion remains kappa>2. The pole is a zero/infinite
spatial-root contribution; it is not a finite-k saddle. Numerically,
direct real-k quadrature gives exp(-t) G00,++(t) close to one: at t=4,
6,8,10 the values are 0.99996321, 1.00002125, 0.99999429, 1.00000113.
8192- and 16384-point grids agree to scaled absolute error below 2e-9
through t=10. Long-time real-k quadrature is not used to establish the
rate because cancellation of growing Fourier components becomes severe.

## 3. Finite OBC: an independent characteristic recurrence

Let p_N(lambda)=det(lambda I-A_N). Eliminating one channel gives

$$
p_N=\det\left[((\lambda+1)I+2S^T)
((\lambda-1)I-2S)+\delta^2I\right].
$$

This tridiagonal determinant has first diagonal entry lambda^2-1+delta^2,
remaining diagonal entries a, upper entries -2(lambda+1) and lower
entries 2(lambda-1). Hence

$$
\begin{aligned}
p_0&=1,\\
p_1&=\lambda^2-1+\delta^2,\\
p_n&=(\lambda^2-5+\delta^2)p_{n-1}
+4(\lambda^2-1)p_{n-2}.
\end{aligned}
$$

We evaluate this recurrence in 70- and 100-digit Decimal arithmetic and
bisect the positive real root for even N. To check that it is the
rightmost root, we independently form the polynomial p_N(s+c) and apply
the Routh root count using 160 and 240 digits. For every reported size
N=4 through 160, the counts to the right of lambda_N-1e-8 and
lambda_N+1e-8 are respectively 1 and 0 at both precisions. This checks
dominance, not just convergence to some characteristic root. These are
precision-refined numerical checks, not interval-arithmetic certificates.

For the even sizes studied,

$$
\Gamma_{\rm OBC}(N)=\lambda_N-\kappa/2,\qquad
\kappa_c(N)=2\lambda_N.
$$

| N | Lossless dominant growth | Critical kappa |
|---|---:|---:|
| 4 | 1.34608857 | 2.69217715 |
| 8 | 1.65020843 | 3.30041685 |
| 16 | 1.89124737 | 3.78249475 |
| 24 | 1.99168672 | 3.98337344 |
| 32 | 2.04662253 | 4.09324506 |
| 48 | 2.10507094 | 4.21014189 |
| 64 | 2.13571102 | 4.27142204 |
| 96 | 2.16735267 | 4.33470533 |
| 160 | 2.19341843 | 4.38683687 |

Strict stability requires kappa above the tabulated threshold. At equality
there is a neutral mode. Odd N need not have a real dominant root; the
recurrence applies to them too, but this particular bisection and table
use even N only.

## 4. Large-chain limit

The recurrence multipliers r satisfy r^2-a r-4(lambda^2-1)=0.
The limiting extended OBC band obeys |r+|=|r-|, equivalently

$$
\lambda^2-5+\delta^2=4i\sqrt{\lambda^2-1}\cos\theta.
$$

Writing y=sqrt(lambda^2-1) gives

$$
y=2i\cos\theta\pm\sqrt{4\sin^2\theta-\delta^2}.
$$

Where the second square root is real, x=lambda^2 has
Re x=5-delta^2-8 cos^2 theta and |x|^2=(5-delta^2)^2-16 cos^2 theta.
Thus Re lambda <= sqrt(5-delta^2), with equality at theta=pi/2.
Where it is imaginary, x is real and at most one, and cannot set the
rightmost edge. The positive real branch above converges to that edge:

$$
\lambda_\infty=\sqrt{5-\delta^2}=2.2338307904,\qquad
\boxed{\kappa_{c,\rm OBC}(\infty)=4.4676615807.}
$$

The recurrence's possible isolated limiting roots arise when a leading
multiplier coefficient vanishes. Substitution of p1 into its multiplier
equation gives 4 delta^2, which is nonzero here, excluding that mechanism.
Additional recurrence calculations give lambda_N=2.20838471, 2.22102824,
2.22740946, 2.23061509 at N=256,512,1024,2048. These track the real branch;
the full Routh counts were performed only through N=160.

The delta->0 and N->infinity limits do not commute. At delta=0 exactly,
p_N=(lambda^2-1)^N and the threshold is 2; arbitrarily small nonzero
coupling can matter after sufficiently long propagation in a large chain.

## 5. SVD topology and its relation to stability

Use Q(k)=omega_ref I-H(k), omega_ref=0, and the Hermitian doubled
matrix [[0,Q],[Q†,0]]. Define nu=(2pi i)^(-1) integral dk partial_k log det Q
with increasing k. The determinant is

$$
\det H(k)=(-2\sin k-i\kappa/2)^2-\delta^2+(1+2\cos k)^2.
$$

For positive kappa its imaginary part vanishes only at k=0,pi. The
point-gap closings occur at kappa=2sqrt(1-delta^2) and
kappa=2sqrt(9-delta^2), yielding

$$
\boxed{\nu(0)=1\quad\text{for}\quad
1.9899748742<\kappa<5.9966657402.}
$$

It is zero outside that interval for kappa>=0 and undefined at its endpoints.
Direct determinant windings on 8192 and 16384 k grids agree with the
formula throughout the loss sweep. Finite-chain singular values are
computed for N=8,16,24,32,48,64 and selected losses. Values near the
double-precision floor are not interpreted as quantitatively accurate
splittings. SVD conditioning and eigenvalue stability are different issues.

The long-chain regimes, away from transition equalities, are:

| Loss interval | Local bulk | Long-chain OBC | PBC | nu(0) |
|---|---|---|---|---|
| 0 to 1.98997487 | Growing | Unstable | Unstable | 0 |
| 1.98997487 to 2 | Growing | Unstable | Unstable | 1 |
| 2 to 4.46766158 | Decaying | Unstable for sufficiently large N | Unstable | 1 |
| 4.46766158 to 5.99666574 | Decaying | Stable | Unstable | 1 |
| Above 5.99666574 | Decaying | Stable | Stable | 0 |

The OBC boundary for a specified finite even N is the tabulated value,
not the limiting one. For example, kappa=4 stabilizes N=16 but not N=32.

The local bulk stability criterion fails to predict finite-chain stability
on 2<kappa<kappa_c(N), entirely within the topological interval. At kappa=3,
the local rate is -0.5, whereas N=16 has rate +0.39124737. Detuning
allows conversion into the oppositely directed channel. Terminating the
chain changes the return paths and the mode spectrum; there is no longer
the exact no-return argument of the decoupled model. No additional
boundary-scattering term was inserted. This demonstrates boundary-sensitive
feedback, without identifying topology with absence of reflection.

## 6. Numerical verification and reproducibility

Run `python src/sweep_quantum2023_detuned.py`. Outputs:

- `data/generated/quantum2023_detuned.json`: parameters, loss grid, winding,
  recurrence and raw eigenvalue comparisons, Routh counts, local-response
  quadrature, singular values, and time-domain traces.
- `figures/quantum2023_detuned.png`: four-panel comparison.
- `notebooks/quantum2023_detuned.ipynb`: loading, plotting and recomputation.

The loss grid is 0 to 8 in steps of 0.05. Exact pi/2 imaginary hopping
avoids an unintended cos(pi/2) term. The Nambu/channel rotation residual
is below 7e-16. Matrix eigenvalues match the recurrence within 1e-7 for
N<=24, but become unreliable as size increases: at N=160 the raw growth
estimate is about 2.51 rather than 2.19341843, despite residuals near 1e-15.
Changing coordinates alone no longer cures this nonnormal sensitivity.

An independent impulse evolution uses scipy expm_multiply on the full
Nambu matrix for N=16, initially particle component j=N/2 equal to one
and all other components zero. Integer times 0 through 500 are sampled.
Fitting log norm on t=420 through 500 gives the rates +0.39124737,
-0.10875263, and -0.60875263 at kappa=3,4,5 within 1e-5. This complex
Nambu impulse is a Green-function column; it need not itself be a physical
conjugate-paired coherent displacement. Physical quadratures span the
same linear response and stability spectrum.

The initial shorter time fit retained oscillatory transients; extending
the time window removed that discrepancy. Fourier quadrature and finite
OBC evolution are deliberately different checks: infinite-chain local
decay does not establish finite-chain stability.

These are new analytical and numerical results for the prescribed
quadratic model, not a reproduction of a mean-field or Gaussian phase diagram.
