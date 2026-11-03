# Rassaert-inspired damped boundaries for the detuned BdG chain

## 1. What is taken from the paper

The published main text of Rassaert et al., *Physical Review Letters* 135,
203603 (2025), Eq. (4), uses a drive envelope epsilon_j=epsilon s_j and
a bare-detuning envelope Delta_j=Delta s_j^2, with

$$
s_j=\tanh\!\left(\frac{d_j}{N_0}\right),\qquad
d_j=\min(j,N+1-j),\quad j=1,\ldots,N.
$$

The loss rate is uniform. This is a smooth suppression of drive and
detuning, not an additional lossy boundary layer. The paper's nonlinear
background changes under this prescription, and so do its effective
Bogoliubov coefficients. Its anomalous coefficients depend on the squared
coherent background; they are not simply the imposed coherent drive.

Sources inspected: the local published PDF `literature/papers/RassaertPRL2025.pdf`,
especially Eqs. (4)-(5), and the accessible supplemental sections III and
VI of [arXiv:2411.08965v1](https://arxiv.org/html/2411.08965v1).
The separate published supplemental link was inaccessible in this session;
the supplemental discussion consulted is explicitly the arXiv version.
[Published article](https://doi.org/10.1103/w32k-x3l5).

## 2. Explicit adaptation to the prescribed quadratic model

We retain the previous bulk model and conventions:
gs=gc=J=1, Delta=0.1, phi=pi/2 exactly, P=0, amplitude damping kappa/2,
and i partial_t Psi=H Psi. Frequencies and rates are in units J=1,
time in 1/J. This is a fixed BdG calculation, with no stationary
mean-field solver and no Gaussian covariance feedback.

There is no independently specified coherent drive epsilon_j in this
quadratic model. Therefore the paper's profile does not uniquely define
a boundary modification of our BdG matrix. Our primary, explicitly
phenomenological adaptation is

$$
\Delta_j=0.1s_j^2,\qquad
g_{s,j}=s_j^2,\qquad
g_{c,j,j+1}=s_js_{j+1}.
$$

This corresponds to the scaling of amplitude-product pairing coefficients
under an assumed amplitude envelope proportional to s_j. We do not infer
that the nonlinear stationary amplitude actually follows the drive, or
that the effective detuning follows the bare detuning. Extending the
profile to the bond pairing gc is our modeling choice; it is not a claim
that this term occurs in the paper's original onsite-only interaction model.

Hopping J, its phase, and kappa remain uniform. Define S[j,j+1]=1,
T=i(S-S^T), D=diag(Delta_j), and real symmetric K with diagonal gs,j
and first off-diagonals gc,j,j+1. Then

$$
H_N=\begin{pmatrix}T+D&K\\-K&-T^*-D\end{pmatrix}
-i\frac{\kappa}{2}I.
$$

Only OBC is used for the tapered chains; no closing bond is inserted.
The homogeneous PBC comparator from the preceding study has the same
directed phase on its closing bond, with total twist -N phi.

A control calculation tapers only Delta_j while retaining gs=gc=1.
The sharp reference has s_j=1 exactly (implemented as width=0).

## 3. Stability calculation and precision checks

Because loss is uniform,

$$
\Gamma_N(\kappa)=\alpha_N-\kappa/2,\qquad
\alpha_N=\max\operatorname{Re}\operatorname{spec}[-iH_N(0)].
$$

Thus the complete loss sweep is obtained by an exact scalar shift; the
critical loss is kappa_c=2 alpha_N. Strict stability requires kappa>kappa_c.
At equality there is no strict exponential decay.

Unstructured double-precision diagonalization is poorly conditioned for
longer chains, as in the preceding homogeneous calculation. We therefore
do not use its apparent convergence as the final stability criterion.
After a constant change of Nambu coordinates and site interleaving, the
lossless generator is real and has bandwidth two. Its onsite block is

$$
A_{jj}=\begin{pmatrix}g_{s,j}&\Delta_j\\-\Delta_j&-g_{s,j}\end{pmatrix},
$$

and the two diagonal intersite blocks are

$$
A_{j,j+1}=\operatorname{diag}(1+g_{c,j},1-g_{c,j}),\qquad
A_{j+1,j}=\operatorname{diag}(g_{c,j}-1,-1-g_{c,j}).
$$

We compute det(lambda I-A) using a banded determinant expansion with
Decimal coefficients. The tanh profiles are also evaluated in Decimal
arithmetic, rather than importing double-precision coefficients. A Routh
count gives the number of roots to the right of a trial real shift.
Bisection locates the spectral abscissa, including complex leading modes;
it does not assume that the dominant eigenvalue is real.

The entire threshold calculation is repeated at 120 and 200 decimal
digits. Both precisions give matching brackets of width below 1e-9 in
growth rate (below 2e-9 in critical kappa). This is a precision-refined
numerical check, not an interval-arithmetic proof. The characteristic
implementation is independently checked against a small dense determinant
and against spectra in the original Nambu coordinates: relative determinant
errors below 6e-15 and spectral discrepancies below 4e-12 at N=8, N0=3.
The sharp result agrees with the earlier independently derived recurrence.

## 4. Results for N=40

| Boundary | N0 | Critical kappa | Growth at kappa=3.5 |
|---|---:|---:|---:|
| Sharp | 0 | 4.16250084 | +0.33125042 |
| Joint pairing/detuning taper | 2 | 3.84542339 | +0.17271170 |
| Joint pairing/detuning taper | 4 | 3.47978130 | -0.01010935 |
| Joint pairing/detuning taper | 6 | 3.07567033 | -0.21216483 |
| Joint pairing/detuning taper | 8 | 2.61912957 | -0.44043521 |
| Detuning-only taper | 8 | 3.72977612 | +0.11488806 |

The joint taper stabilizes this chain at kappa=3.5 for N0=4,6,8.
The N0=8 joint taper also stabilizes kappa=3, with growth -0.19043521,
whereas the sharp N=40 chain grows at +0.58125042.

The original homogeneous local bulk criterion is kappa>2; the original
PBC threshold is 5.99666574. These remain reference values for the
original uniform operator, not exact criteria for the tapered chain.
Even N0=8 leaves a threshold above 2: the profile helps, but does not
restore the homogeneous local criterion exactly.

## 5. How much of the bulk is changed?

A tanh profile has tails, so at finite N it never equals its asymptotic
bulk value exactly. This matters when interpreting a very wide taper.

| N=40, N0 | Central gs (bulk value 1) | Central Delta | Sites with s_j^2>0.99 |
|---:|---:|---:|---:|
| 2 | 0.99999999 | 0.10000000 | 30 |
| 4 | 0.99981842 | 0.09998184 | 18 |
| 6 | 0.99492240 | 0.09949224 | 6 |
| 8 | 0.97340777 | 0.09734078 | 0 |
| 12 | 0.86696510 | 0.08669651 | 0 |
| 16 | 0.71958513 | 0.07195851 | 0 |

The N0=4 example is therefore the cleanest stabilization result with
nearly unchanged central coefficients: their fractional reduction is
only 0.0182%. For N0=8 the central reduction is already 2.66%.
Thresholds below 2 for very broad tapers should not be interpreted as
beating the stability bound of an unchanged infinite bulk: most or all
of those finite chains have weaker pairing than the original model.

## 6. Chain-length dependence

High-precision critical loss for the joint taper:

| N | N0=2 | N0=4 | N0=6 | N0=8 |
|---:|---:|---:|---:|---:|
| 16 | 3.007415 | 1.884525 | 1.080220 | 0.806300 |
| 32 | 3.700244 | 3.234817 | 2.690818 | 2.106181 |
| 40 | 3.845423 | 3.479781 | 3.075670 | 2.619130 |
| 64 | 4.070199 | 3.843144 | 3.609556 | 3.361557 |
| 80 | 4.147244 | 3.965159 | 3.780781 | 3.589276 |

At kappa=3.5, N0=8 stabilizes N=64 but not N=80. Increasing N0 to 12
at N=80 gives kappa_c=3.17520833 and stabilizes it; its central gs is
approximately 0.9949224. A fixed-width boundary does not provide a
size-independent stability guarantee. No extrapolated threshold at
infinite N is claimed for these tapered profiles.

The detuning-only control also lowers thresholds, but more weakly. At
N=40 and N0=4,6,8 its critical losses are 3.81486873, 3.76404899,
3.72977612, respectively. All remain unstable at kappa=3.5. Reducing
pairing near the edges is consequential for the stabilization reported here.

## 7. Independent time evolution and SVD

At N=40, kappa=3.5 we evolve a physical coherent perturbation with equal
particle and hole amplitudes 1/sqrt(2) at zero-based site N/2, all others
zero. scipy expm_multiply evolves the full Nambu generator through Jt=2000
with spacing 2. A known scalar growth factor is removed during evolution
and restored in log norm to avoid overflow and underflow. This is an
exact change of normalization, not a projection onto a selected mode.

Fits to log norm over Jt=1000 to 2000 give:

| Boundary | N0 | Spectral prediction | Time-domain fit |
|---|---:|---:|---:|
| Sharp | 0 | +0.33125042 | +0.33125136 |
| Joint taper | 4 | -0.01010935 | -0.01010992 |
| Joint taper | 6 | -0.21216483 | -0.21216601 |
| Joint taper | 8 | -0.44043521 | -0.44043311 |
| Detuning only | 8 | +0.11488806 | +0.11492307 |

Differences are below 4e-5. Earlier time windows have larger deviations
from oscillatory transients. The script allows 2e-3 for this independent
fit; it does not use time-fit precision to set the stability boundary.

The data also record the two lowest singular values of Q=-H at reference
frequency zero. They remain strongly separated in the examples, including
both stable and unstable ones. For example, the stable joint N0=6 case
has smin about 1.6e-9 and the next singular value about 0.816. This does
not make dynamical stability a consequence of SVD topology. The original
homogeneous bulk winding is +1 at kappa=3.5; a translation-invariant
winding is not assigned to the entire inhomogeneous chain. No untested
slowly-varying local-invariant approximation is used to certify stability.
Singular values near the double-precision floor are qualitative diagnostics.

## 8. Interpretation and reproducibility

This prescription can stabilize the Delta=0.1 BdG chain without adding
loss. It modifies the boundary-region amplification and conversion
processes. The calculation does not extract a reflection coefficient,
so it does not establish that the taper is reflectionless, nor isolate
reduced reflection from reduced boundary gain as the mechanism.

This is a new calculation for an explicitly specified quadratic adaptation.
It is not a reproduction of the paper's nonlinear Gaussian phase diagram,
and does not determine which mean-field solution of a microscopic
cross-Kerr model would generate these coefficients.

Reproduction: `python src/bdg_damped_boundaries.py` from the project root.
Reusable precision routines are in `src/banded_stability_decimal.py`.
The initial double-precision survey uses N=16,24,32,40,64,80, widths
0,1,2,4,6,8 and both profile types. The refined survey uses N=16,32,40,64,80,
widths 2,4,6,8 and both types, plus joint widths 12,16 at N=40,64,80.
Raw results are preserved as conditioning diagnostics, not used for final
thresholds. There are no random parameters or stationary initial guesses.

Outputs: `data/generated/bdg_damped_boundaries.json`,
`figures/bdg_damped_boundaries.png`, and the loading/recomputation notebook
`notebooks/bdg_damped_boundaries.ipynb`.

The extension to N=128,256,512, with fixed and size-scaled boundary widths,
is documented in `bdg_long_boundary_scaling.md` and
`notebooks/bdg_long_boundaries.ipynb`.
