# N=100: sharp left edge and a right-only absorbing region

## Question and result

Can the original homogeneous bulk threshold kappa_c=2 be recovered by
changing only the right boundary of the Delta=0.1 chain? In the finite
families tested here, **no**. Smooth profiles can lower the finite-chain
threshold substantially, but none makes it coincide with 2. This is a
negative result for a documented search, not a proof that no engineered
absorbing boundary can work.

The distinction between a low reflection amplitude at selected frequencies
and global dynamical stability is essential in this strongly amplifying
system. We calculate both independently below.

## 1. Fixed model and exactly unchanged left region

N=100, bulk gs=gc=J=1, Delta=delta=0.1, phi=pi/2, collective pump P=0,
hbar=1, lattice spacing one. Evolution is i partial_t Psi=H Psi, with
uniform baseline amplitude damping kappa/2. Positive Im omega means growth.
There is a sharp OBC at the left end and no closing bond. These are fixed
quadratic coefficients; no nonlinear mean-field or Gaussian state is solved.

Let L be the number of modified right sites and M=100-L the unchanged
core length. For j<=M, all onsite parameters, loss, and bonds entirely
within the core remain exactly equal to the original model. For j>M,
put x_j=(j-M)/L. Unlike an untruncated tanh tail, this construction does
not weaken the left edge or the retained core.

For the joint pairing/detuning taper,

$$
\Delta_j=0.1s_j^2,\qquad g_{s,j}=s_j^2,\qquad
g_{c,j,j+1}=s_js_{j+1},
$$

where s=1 in the core. The five initial right-layer shapes are

$$
\begin{aligned}
s_{\rm linear}(x)&=1-x,\\
s_{\rm cubic}(x)&=1-3x^2+2x^3,\\
s_{\rm quintic}(x)&=1-10x^3+15x^4-6x^5,\\
s_{\tanh,a}(x)&=\tanh[a(1-x)]/\tanh(a),\qquad a=2,4.
\end{aligned}
$$

The quintic profile joins the core with two continuous zero derivatives.
At the last site s=0. Hopping J and its phase remain unchanged. These
are right-only extensions of our previous prescribed quadratic adaptation
of Rassaert's boundary smoothing, not profiles derived from the paper's
nonlinear stationary solution.

We also test added-loss layers, separately from the above prescription.
They leave gs=gc=1 and Delta=0.1 everywhere, and introduce extra amplitude
damping beta_j in the right region. The total local Lindblad rate is
kappa_j=kappa+2 beta_j. We test beta=eta x^p for p=1,2,4,8, as well as
beta=eta[1-s_quintic(x)]. Combined taper/loss cases use the quintic s and
the latter beta. No additional loss acts in the left core.

## 2. Stability search

In real rotated and site-interleaved coordinates, the generator A=-iH
has onsite and nearest-neighbour blocks

$$
A_{jj}=\begin{pmatrix}g_{s,j}-\beta_j&\Delta_j\\-\Delta_j&-g_{s,j}-\beta_j\end{pmatrix},
$$

$$
A_{j,j+1}=\operatorname{diag}(1+g_{c,j},1-g_{c,j}),\qquad
A_{j+1,j}=\operatorname{diag}(g_{c,j}-1,-1-g_{c,j}).
$$

Baseline loss is then subtracted as kappa I/2. Consequently

$$
\Gamma_{100}(\kappa)=\alpha-\kappa/2,\qquad
\kappa_c=2\alpha,\qquad \alpha=\max\operatorname{Re}\operatorname{spec}A.
$$

The original local-bulk rate is 1-kappa/2. All thresholds below are
obtained from high-precision characteristic polynomials and full Routh
root counts, not from unstructured double-precision eigenvalues. Working
precision is 160 decimal digits, with bracket endpoints checked again
at 260 digits. Critical-loss brackets have width below 2e-5; positive
counts below and zero counts above verify spectral ordering, including
cases with complex dominant modes. These are precision-refined numerical
checks rather than interval-arithmetic certificates.

There are 73 distinct cases: a sharp reference; five taper families at
L=20,40,60,80,90; quintic extensions at L=92,94,96,98,99; smoothstep
loss-only layers at L=40,80 with eta=1,3,10; polynomial loss-only layers
at L=80,90 with p=1,2,4,8 and eta=1,3,10; and combined quintic/loss
layers at L=60,80,90 with eta=0.3,1,3,10. This is a finite family/grid
search, not global optimization over all spatial coefficients.

Representative thresholds for the quintic taper, with no added loss:

| Right-layer L | Unchanged sites M | Critical kappa |
|---:|---:|---:|
| 0 | 100 | 4.33985 |
| 20 | 80 | 4.16174 |
| 40 | 60 | 4.01417 |
| 60 | 40 | 3.77204 |
| 80 | 20 | 3.17872 |
| 90 | 10 | 2.52729 |
| 96 | 4 | 2.21226 |
| 98 | 2 | 2.14884 |
| 99 | 1 | 2.11462 |

With at least 20 unchanged sites, the lowest threshold among all tested
cases is 3.04964, for L=80 and beta=10x^4 (loss only). With at least
10 unchanged sites, the best is the quintic L=90 taper at 2.52729.
The loss-only L=90, beta=10x^2 case reaches 2.67590. Adding extra loss
to the quintic taper did not lower its threshold on the tested grid.
Larger loss is therefore not automatically better.

At kappa=2.1, **all 73 tested chains are unstable**, despite decay of
the original infinite-chain local response. Reducing the threshold by
tapering almost all sites does not demonstrate a boundary treatment that
preserves a substantial original bulk region.

## 3. Modal reflection, independently of the stability calculation

We measure returning/incident amplitude at the entrance of the right
layer using the two spatial modes of the unchanged bulk. The calculation
uses unit Euclidean norm of each two-component mode vector. This is not
a flux-normalized power reflection coefficient; gain permits magnitudes
above one. For an unstable terminated system, the frequency-domain
quantity is a formal scattering diagnostic, not an asserted stationary
response to a drive switched on forever.

For a real frequency omega set u=kappa/2-i omega. The lossless bulk
generator in these real channel coordinates is

$$
A(z)=\begin{pmatrix}1+2z&\delta\\-\delta&-1-2/z\end{pmatrix},\quad \delta=0.1.
$$

Its spatial polynomial is

$$
-2(u+1)z^2+(u^2-4.99)z+2(u-1)=0.
$$

The outgoing right-going root is selected by continuation from large
positive Re u, where z_R approximately equals -2/u. In Re u>0 use

$$
F^{1/2}(u)=\sqrt{u^2+2.61}\sqrt{u^2+3.41},\qquad
z_R=\frac{u^2-4.99-F^{1/2}}{4(u+1)},\quad
z_L=-\frac{u-1}{(u+1)z_R}.
$$

The square roots have positive real part on this half-plane. Spatial
roots are labeled by this causal continuation, not simply by whether
their modulus is smaller than one at the final loss.

Corresponding unnormalized vectors can be taken as

$$
v_R=\begin{pmatrix}\delta/(u-1-2z_R)\\1\end{pmatrix},\qquad
v_L=\begin{pmatrix}1\\-\delta z_L/[(u+1)z_L+2]\end{pmatrix}.
$$

After normalizing each vector, write the uniform-core solution as a
superposition with incident amplitude one and returning amplitude R.
The amplitudes are referenced by extrapolation to the last uniform site.
Include that site and all L modified sites in the boundary solve.
On the bond from the preceding unchanged site, only p propagates left
and q propagates right. Eliminating the external fields gives

$$
p_0=v_{R,p}+R v_{L,p},\quad
q_{-1}=v_{R,q}/z_R+R v_{L,q}/z_L,
$$

$$
\chi=\frac{v_{L,q}}{z_Lv_{L,p}}
=-\frac{\delta}{(u+1)z_L+2}.
$$

For B=uI-A_layer, add 2chi to its (q0,p0) entry and put
-2(v_R,q/z_R-v_R,p chi) in the q0 source entry. Solve this banded
system, then extract R=(p0-v_R,p)/v_L,p. This includes the modified
interface bond and all propagation/conversion inside the layer.

Banded complex solves are repeated at 100 and 160 decimal digits.
The reported double representations agree to relative/absolute tolerance
1e-12 using max(1,|R|) as the scale. Naive double-precision solves fail
badly for some near-zero-frequency responses and are not used for the
final reflection values.

At kappa=2.1, representative returning amplitudes are:

| Right termination | abs R at omega=0 | abs R at omega=0.5 | abs R at omega=1 |
|---|---:|---:|---:|
| Sharp | 1.9427 | 0.19856 | 0.10612 |
| Quintic L=80 | 24.2428 | 3.787e-6 | 7.793e-6 |
| Quintic L=90 | 48.4765 | 2.739e-6 | 5.460e-6 |
| Quintic L=96 | 73.3444 | 2.465e-6 | 4.493e-6 |
| Loss only, L=90, beta=10x^2 | 0.03508 | 4.068e-5 | 8.804e-5 |

Thus the smooth taper is highly effective at some sampled frequencies,
but emphatically not a nearly perfect broadband absorber. The narrower
frequency window near zero is crucial here. The saved frequency grid is
0,0.001,0.003,0.01,0.03,0.05,0.1,0.2,0.3,0.5,0.75,1,1.5,2,3.
Interpolated plot lines are guides between samples, not a continuous-band
bound on reflection.

## 4. Why a tiny return can still destabilize the chain

The sharp left boundary imposes q_0=0. With M unchanged sites, the
uniform-core transfer therefore requires

$$
R_{\rm required}(u)=-\frac{v_{R,q}}{v_{L,q}}
\left(\frac{z_L}{z_R}\right)^M.
$$

The eigenmode/feedback condition is R(u)=R_required(u), including phase.
The associated round-trip magnitude is |R/R_required|. A magnitude
greater than one alone is not a stability proof; the independent
characteristic-root calculation establishes instability. Possible poles
of subsystem responses must also be considered in a full feedback argument.

For M=10, kappa=2.1 and omega=0,

$$
|R_{\rm required}|\simeq1.048\times10^{-13}.
$$

The quintic L=90 value is approximately 48.48, and even the loss-only
L=90 value 0.03508 is vastly larger than this scale. For longer unchanged
cores the required return suppression near the bulk threshold becomes
still more stringent. This explains why visually small reflection or
small reflection at omega=0.5 is not enough.

In the present hopping convention the most amplifying Bloch modes near
k=0 have real group velocity -2J: their amplification travels toward
the sharp LEFT end. A right absorber must suppress the subsequent return
path extremely well; it does not remove that first sharp left boundary.

## 5. An exact transparent benchmark exists, but is not a static profile

Replacing the right end by an outgoing semi-infinite continuation amounts
to a frequency-dependent generator self-energy

$$
\Sigma_R(u)=2z_R(u)
$$

on the last p-coordinate diagonal. It enforces p_(N+1)=z_R p_N.
For a uniform reference segment this gives R=0 by construction. Our
implementation reproduces |R|<1e-85 at the four test frequencies
0,0.5,1,2 with 100-digit arithmetic. A second check reproduces the sharp
N=8 boundary-mode condition at its independently known critical loss
3.300416850278589 to relative error below 1e-10.

This identifies a concrete ideal load against which to judge candidate
absorbers. Its frequency dependence corresponds to memory in the time
domain. It is **not** a frequency-independent local damping profile, nor
is it one of the finite 200-by-200 BdG matrices tested above. We have not
implemented a dynamical auxiliary-mode approximation to this load here,
and do not report its zero reflection as a successful finite-profile result.

## Files and reproduction

- `src/bdg_right_absorber.py`: resumable 73-case spectral search.
- `src/bdg_right_reflection.py`: double and precision-refined modal matching;
  only the Decimal results are used for final reflection claims.
- `src/decimal_complex_banded.py`: complex Decimal banded solver.
- `src/analyze_right_absorber.py`: reflection checks and comparison figure.
- `data/generated/bdg_right_absorber.json`, `bdg_right_reflection.json`:
  numerical parameters, root-count brackets, and reflection diagnostics.
- `figures/bdg_right_absorber.png` and
  `notebooks/bdg_right_absorber.ipynb`: presentation and recomputation.

Run the spectral search first and `python src/analyze_right_absorber.py`
second. Cached spectral cases are retained; remove or rename the dataset
explicitly to force a full rerun. No random initial conditions are used.
These results are for the same fixed quadratic model with documented
right-boundary modifications, not a nonlinear stationary-state calculation.
