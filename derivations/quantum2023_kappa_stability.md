# Loss sweep in the Quantum 2023 BdG convention

## 1. Source and conventions

Source: A. Gómez-León, T. Ramos, A. González-Tudela and D. Porras,
*Driven-dissipative topological phases in parametric resonator arrays*,
[Quantum 7, 1016 (2023)](https://doi.org/10.22331/q-2023-05-23-1016),
[arXiv version](https://arxiv.org/pdf/2207.13715).
We use the **dynamical-matrix definitions**, Eqs. (7)-(10) and (28)-(31).
The paper calls the local Lindblad loss rate gamma; here it is kappa.
Collective pumping P is set to zero. Units are hbar=J=1, lattice spacing one.
These are fixed quadratic parameters, not a mean-field stationary-state calculation.

The onsite anomalous BdG element is gs and the anomalous element on each bond
is gc. The Hamiltonian onsite term implementing this matrix is
gs (a_j^2 + a_j^{dagger 2})/2 for real gs. The displayed Eq. (2) in the
arXiv paper writes gs without this half, whereas Eq. (10) uses gs as the
BdG element. We follow Eq. (10) and the explicit Bloch matrix, rather than
silently introducing a factor two. The hopping orientation is also fixed by
Eq. (9) and the Bloch matrix below.

For Fourier convention u_j proportional to exp(ikj), and evolution
i partial_t Psi = H Psi, the matrix is

$$
H(k)=\left[-2J\sin\phi\sin k-i\frac{\kappa}{2}\right]I
+\left[\Delta+2J\cos\phi\cos k\right]\sigma_z
+i[g_s+2g_c\cos k]\sigma_y.
$$

Stability means max Im omega < 0, because a mode evolves as exp(-i omega t).
We study gs=gc=J=1, Delta=0, phi=pi/2 exactly. For PBC the same complex
hopping is used on the closing bond; the directed total phase is -N phi
in the lower-diagonal convention of Eq. (9), modulo 2pi. No boundary
drive or detuning taper is introduced.

## 2. A constant rotation separates two one-way channels

Since the sigma_z term vanishes identically, use normalized eigenvectors
(1,i)/sqrt(2) and (1,-i)/sqrt(2) of sigma_y. This is a unitary change of
coordinates for analyzing the linear equation, not a claim about a canonical
bosonic Bogoliubov transformation. With z=exp(ik), the eigenvalues are

$$
\begin{aligned}
\omega_+(k)&=-2\sin k-i\kappa/2+i(1+2\cos k)
=i(1-\kappa/2)+2iz,\\
\omega_-(k)&=-2\sin k-i\kappa/2-i(1+2\cos k)
=-i(1+\kappa/2)-2iz^{-1}.
\end{aligned}
$$

Let S have entries S[j,j+1]=1, j=0,...,N-2. Thus (Su)_j=u_{j+1}.
For OBC, S^N=0. The two generators A=-iH are

$$
A_+=(1-\kappa/2)I+2S,\qquad
A_-=-(1+\kappa/2)I-2S^T.
$$

Each channel propagates in one direction, and the two channels do not mix,
including at the open ends. This is the reason boundary feedback is absent
in this particular clean, fine-tuned example.

## 3. Convective propagation, without a finite saddle approximation

Put a unit perturbation in channel + at site m. The propagator to n<=m is

$$
G^+_{nm}(t)=e^{(1-\kappa/2)t}\frac{(2t)^{m-n}}{(m-n)!},
$$

and it vanishes for n>m. For channel -,

$$
G^-_{nm}(t)=e^{-(1+\kappa/2)t}\frac{(-2t)^{n-m}}{(n-m)!},\quad n\geq m,
$$

and it vanishes for n<m. In the original Nambu coordinates,

$$
G^{\mathrm{OBC}}_{nm}(t)=U
\begin{pmatrix}G^+_{nm}(t)&0\\0&G^-_{nm}(t)\end{pmatrix}U^\dagger,
\quad U=\frac{1}{\sqrt2}\begin{pmatrix}1&1\\i&-i\end{pmatrix}.
$$

At any fixed separation, the polynomial factor does not change the
exponential rate. The maximum fixed-position rate is therefore

$$
\Gamma_{\mathrm{bulk}}(v=0)=1-\kappa/2.
$$

The same coefficients follow from extracting the Laurent coefficients of
exp[-i H(k)t] on the unit circle. Here the one-way limit sends the usual
zero-velocity saddle to z=0 or infinity: differentiating omega_+ gives
-2z and differentiating omega_- gives -2/z. There is no finite saddle to
insert into a generic saddle-point formula. The coefficient calculation
is exact and avoids this singular limit.

To see convection explicitly, let the distance along channel + be d=vt.
Stirling's formula yields

$$
\Gamma_+(v)=1-\kappa/2+v[1+\log(2/v)],\qquad v>0.
$$

Its maximum is at v=2, with growth rate 3-kappa/2. Consequently, for
2<kappa<6 a localized disturbance decays at a fixed position while it can
grow along a moving ray. This is convective amplification with local decay.
The minus-channel moving-ray maximum is 1-kappa/2.

## 4. Exact OBC and PBC stability

At fixed finite N the Taylor series terminates:

$$
e^{A_+t}=e^{(1-\kappa/2)t}\sum_{d=0}^{N-1}\frac{(2t)^dS^d}{d!}.
$$

Thus every initial perturbation decays for kappa>2. Equivalently, the
OBC eigenvalues are i(1-kappa/2) and -i(1+kappa/2), each with algebraic
multiplicity N. Except for N=1 these are nontrivial Jordan blocks.

$$
\boxed{\Gamma_{\mathrm{OBC}}=1-\kappa/2,\quad
\text{strict finite-OBC stability iff }\kappa>2.}
$$

At kappa=2 there is polynomial growth for N>1: zero exponential rate
does not imply bounded dynamics. For N=1 the plus channel is neutral.
There can be very large transients even above threshold. The claim is
fixed-finite-N, long-time stability, not a bound uniform in system size.

PBC makes S cyclic, allowing repeated circulation. The largest imaginary
part of the Bloch spectrum occurs in channel + at k=0, available for every N:

$$
\boxed{\Gamma_{\mathrm{PBC}}=3-\kappa/2,\quad
\text{strict PBC stability iff }\kappa>6.}
$$

At kappa=6 PBC has a neutral mode. Open-chain boundaries cause no reflected
return path here: both directed channels terminate without converting into
one another. Changing Delta, phi, couplings, or boundary mixing invalidates
this special reduction and requires a new stability calculation.

## 5. SVD topology at a specified reference frequency

Use Q(k)=omega0 I-H(k) and the Hermitian doubled matrix [[0,Q],[Q†,0]].
For real omega0 define nu=(2pi i)^(-1) integral dk partial_k log det Q.
Increasing k from -pi to pi, the plus circle has positive orientation and
the minus circle negative orientation. Each has radius 2, giving

$$
\nu(\omega_0)=
\mathbf1\{\omega_0^2+(\kappa/2-1)^2<4\}
-\mathbf1\{\omega_0^2+(\kappa/2+1)^2<4\}.
$$

On either equality the point gap closes and the invariant is undefined.
At omega0=0 and nonnegative kappa:

| Loss | Zero-frequency winding | OBC | PBC |
|---|---|---|---|
| 0 <= kappa < 2 | 0 (two opposite channel windings) | Unstable | Unstable |
| kappa = 2 | Undefined, bulk point gap closes | Polynomial growth, N>1 | Unstable |
| 2 < kappa < 6 | +1 | Stable, possibly large transients | Unstable |
| kappa = 6 | Undefined, bulk point gap closes | Stable | Neutral mode |
| kappa > 6 | 0 | Stable | Stable |

For 0<=kappa<2, two exponentially small singular values can occur despite
zero total winding: the opposite contributions cancel. They are not
protected by the full determinant invariant against channel mixing.
For 2<kappa<6 there is one protected small singular value. Writing
a=kappa/2-1 in (0,2), its large-N asymptotic is

$$
s_{\min}\sim 2\left[1-(a/2)^2\right](a/2)^N.
$$

This follows from the nearly null geometric vector of the bidiagonal Q+
block and its oppositely localized left vector. At kappa=6 decay is
algebraic (bulk gap closes); above 6 the minimum tends to kappa/2-3.

The bulk criterion is exact throughout this parameter line, in both the
trivial and nontrivial regions. Its validity is therefore **not caused by
nonzero winding**. The exact decoupling proves it. This does not reverse
the earlier phase-sweep counterexample with nonzero detuning.

## 6. Numerical comparison and precision

Reproduce with `python src/sweep_quantum2023_kappa.py`; the notebook
`notebooks/quantum2023_kappa.ipynb` loads and plots the saved results.
Data: `data/generated/quantum2023_kappa.json`; figure:
`figures/quantum2023_kappa.png`. No random initialization or stationary
solver is used. Loss grid: 0 to 8 in steps of 0.1. Winding grids: 4096
and 8192 points, excluding analytically identified gap closings.

Numerical checks include channel spectra, full Nambu-to-channel matrix
rotation, PBC spectra for N=16,32,64, finite SVD for N=8,16,32,64,96,
and full-matrix expm versus the finite polynomial propagator at N=8.
The rotation residual is below 9e-16 in infinity norm. For t<=10 the
relative propagator discrepancy is below 2.2e-15. Longer-time discrepancies
are retained in the dataset rather than used to infer a new stability
threshold. Acceptance tolerances are 1e-13 (rotation), 1e-12 (spectral
identities), and 1e-7 (short-time full-basis propagator comparison).

At kappa=3 the exact/channel-preserving rate is -0.5 for every N, but
unstructured double-precision diagonalization of the full Nambu matrix gives:

| N | Raw Nambu eigenvalue estimate | Exact rate |
|---|---|---|
| 8 | -0.48385 | -0.5 |
| 16 | -0.33049 | -0.5 |
| 32 | +0.08865 | -0.5 |
| 64 | +0.57454 | -0.5 |
| 96 | +0.83063 | -0.5 |

These are spurious instabilities of an exactly defective matrix. Small
backward errors can split an N-fold Jordan eigenvalue by a scale proportional
to error^(1/N). A tiny residual alone does not certify accurate eigenvalues.
The results depend on numerical library/precision and are diagnostic, not
physical observations. Exact pi/2 hopping is represented by imaginary
entries, avoiding a spurious cos(pi/2) term; this alone does not cure an
unstructured eigensolver. The constant channel rotation and nilpotent
propagator supply independent algebraic checks.

The paper discusses size-dependent stability in Sec. 3.2.1 / Fig. 6. We do
not claim to reproduce that boundary: the exact clean matrix specified here
has the size-independent spectrum proved above. The observed numerical
sensitivity is a concrete reason to investigate precision and perturbations
before interpreting an apparent size dependence at this endpoint. We have
not audited the authors' code or established the origin of their figure.

This derivation and the precision comparison are new calculations using
the paper's BdG conventions, not a reproduction of its complete phase diagram.
