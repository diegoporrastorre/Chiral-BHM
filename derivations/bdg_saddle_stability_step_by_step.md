# Bulk convective stability, open boundaries, and SVD topology

This document derives a fixed-background stability diagnostic for an onsite-paired dissipative bosonic chain and explains the additional assumptions needed to use it for OBC. The accompanying phase sweep is reported in `docs/bdg_phi_stability_results.md`. All display mathematics uses `$$` delimiters.

## 1. Model and conventions

Take lattice spacing and hbar equal to one. Sites are j=1,...,N. The prescribed homogeneous quadratic Hamiltonian is

$$
H_q=\sum_j\left[\tilde\Delta b_j^\dagger b_j+
\frac12(g_s b_j^{\dagger 2}+g_s^*b_j^2)\right]
+J\sum_{j=1}^{N-1}\left(e^{i\phi}b_{j+1}^\dagger b_j+e^{-i\phi}b_j^\dagger b_{j+1}\right).
$$

Here J>0, tilde Delta and phi are real, and g_c=0 means no intersite anomalous pairing. A uniform phase rotation makes g_s=g>=0 without changing phi. Local loss is kappa D[b_j], with D[b]rho=b rho b^dagger-{b^dagger b,rho}/2. The amplitude-loss rate is gamma=kappa/2.

OBC means hard truncation; there is no damping profile. PBC adds J exp(i phi) b_1^dagger b_N+h.c., giving total twist N phi modulo 2pi. Parameters are held fixed when changing boundary conditions.

In the Nambu ordering (b_1,...,b_N,b_1^dagger,...,b_N^dagger),

$$
i\partial_t\Psi=\mathcal H\Psi,\qquad
\mathcal H=\begin{pmatrix}A&G\\-G^*&-A^*\end{pmatrix}-i\gamma I,
\quad G=g_s I,
$$

with A_jj=tilde Delta, A_(j+1,j)=J exp(i phi). Input noise does not change the deterministic stability eigenvalues. Since a mode evolves as exp(-i omega t), its growth rate is Im(omega). Strict finite-system stability requires every eigenvalue to have negative imaginary part. At zero growth, Jordan blocks must also be checked.

This is not a self-consistent Gaussian calculation. For the V=0 mean-field expansion of the paper's U n_j^2 convention, the coefficients would be

$$
\tilde\Delta_j=\Delta_j+U+4U|\alpha_j^\star|^2,
\qquad g_{s,j}=2U(\alpha_j^\star)^2.
$$

We instead prescribe homogeneous coefficients to isolate the spectral effect of boundaries. No stationary mean-field branch or residual is asserted for these examples.

## 2. Bloch matrix and temporal (PBC) instability

With b_j proportional to exp(ikj), and spinor (b_k,b_-k^dagger), define

$$
a(k)=2J\sin\phi\sin k,\qquad
d(k)=\tilde\Delta+2J\cos\phi\cos k.
$$

Then

$$
\mathcal H(k)=[a(k)-i\gamma]I+B(k),\quad
B(k)=\begin{pmatrix}d(k)&g_s\\-g_s^*&-d(k)\end{pmatrix},
\quad B(k)^2=[d(k)^2-|g_s|^2]I.
$$

Thus

$$
\omega_\pm(k)=a(k)-i\gamma\pm\sqrt{d(k)^2-|g_s|^2}.
$$

For real k, a and d are real. The continuous-momentum maximum growth is

$$
\Gamma_{\rm PBC}=-\gamma+\sqrt{\max(0,|g_s|^2-d_{\min}^2)},
\quad d_{\min}=\max(0,|\tilde\Delta|-2|J\cos\phi|).
$$

A finite ring samples k=2pi l/N; its maximum can differ. The continuous result is exact for the bulk, not automatically for every finite ring.

## 3. Start from a localized perturbation

The infinite-chain site propagator is a 2x2 matrix:

$$
\mathcal G^\infty_{nm}(t)=\int_{-\pi}^{\pi}\frac{dk}{2\pi}
e^{ik(n-m)}e^{-i\mathcal H(k)t}.
$$

It maps initial particle and hole amplitudes at site m into those at n. Away from r(k)=0, where r(k)=sqrt(d(k)^2-|g_s|^2), the spectral projectors are

$$
P_\pm(k)=\frac12\left[I\pm\frac{B(k)}{r(k)}\right],
\quad e^{-i\mathcal H(k)t}=\sum_{\sigma=\pm}P_\sigma(k)e^{-i\omega_\sigma(k)t}.
$$

For numerical evaluation, avoid the singular-looking projectors by using the exact entire matrix expression

$$
e^{-i\mathcal H(k)t}=e^{-ia(k)t-\gamma t}
\left[\cos(r(k)t)I-i\frac{\sin(r(k)t)}{r(k)}B(k)\right],
$$

where sin(rt)/r tends to t at r=0. A unit excitation in each Nambu component gives the two columns of the impulse response. Physical classical perturbations impose the usual conjugacy between components; the full matrix norm supplies an upper bound for them.

## 4. Longhi's moving-ray argument

Set v=(n-m)/t. The exponential in each band integral is exp[it(kv-omega_sigma(k))]. Its stationary points satisfy

$$
\omega_\sigma'(k_s)=v.
$$

k_s is generally complex. For an accessible, nondegenerate saddle, steepest descent gives a contribution of the form

$$
\mathcal G^\infty_{nm}(t)\sim t^{-1/2}C_\sigma(v)
e^{it[k_s v-\omega_\sigma(k_s)]}.
$$

The exponential rate is Im[omega_sigma(k_s)-v k_s]. Multiple equal-rate saddles can interfere. Degenerate saddles change the algebraic prefactor. The exponent is best defined by a limsup of log(norm)/t to avoid zeros caused by oscillation.

At fixed n-m, v tends to zero. Write omega_sigma=Omega_sigma-i gamma and define

$$
\Gamma_{\rm abs}=\max_{\text{contributing saddles}}
\operatorname{Im}\Omega_\sigma(k_s),\qquad \Omega_\sigma'(k_s)=0.
$$

Uniform loss does not change the saddle locations, so the fixed-position rate is Gamma_abs-gamma. Accordingly,

$$
\boxed{\gamma>\Gamma_{\rm abs}\quad\text{or}\quad\kappa>2\Gamma_{\rm abs}}
$$

ensures exponential local bulk decay. If Gamma_PBC>0 at the same time, the bulk is convectively unstable: disturbances grow on some moving rays while decaying at fixed sites. This is not decay of the total norm on an infinite chain. The equality case needs the algebraic behavior and possible degeneracies.

## 5. Compute stationary saddle candidates and select them

Let s=J sin(phi), c=J cos(phi). The unsquared equation is

$$
2s\cos k_s\mp\frac{2c\,d(k_s)\sin k_s}{\sqrt{d(k_s)^2-|g_s|^2}}=0.
$$

Set z=exp(ik_s), Q(z)=tilde Delta z+c(z^2+1). Squaring and clearing denominators yields

$$
s^2(z^2+1)^2[Q(z)^2-|g_s|^2z^2]
+c^2Q(z)^2(z^2-1)^2=0.
$$

Check nonzero roots against the original derivative to remove artifacts. Square-root degeneracies need separate treatment.

Not all surviving algebraic saddles contribute. The spatial dispersion polynomial at zero uniform loss is

$$
F(z,\Omega)=[\Omega z+is(z^2-1)]^2-Q(z)^2+|g_s|^2z^2.
$$

At a saddle, two spatial roots coalesce. The sweep follows that pair along Omega=Omega_s+i y to large positive y, where the four roots divide into two small-|z| and two large-|z| roots. Coalescence between roots coming from opposite sets is the numerical pinching diagnostic. The initial separation and tracking resolution are refined independently.

This continuation is a diagnostic along a specified path, not a universal implementation of contour topology. Other branch points on the path, root-assignment ambiguities, or degeneracies require separate analysis. Compare with direct Fourier propagation and exact limiting cases. Never use the maximum imaginary part of every polynomial root indiscriminately: it can exceed the maximum real-k growth, contradicting a bound on the original integral.

## 6. What this does and does not establish for OBC

The exact finite-chain propagator is

$$
\mathcal G^{\rm OBC}(t)=e^{-i\mathcal H_{\rm OBC}t}
=\frac{1}{2\pi i}\oint_\Gamma e^{-i\omega t}
(\omega I-\mathcal H_{\rm OBC})^{-1}d\omega.
$$

The contour surrounds the spectrum counterclockwise. This includes all boundary scattering. The bulk criterion predicts OBC stability only if no boundary mode or reflected feedback produces a larger growth rate. This must be tested rather than inferred from drift alone. Finite-N eigenvalues are compared with the bulk threshold in the sweep, and conditioning is recorded because small residuals do not ensure accurate nonnormal eigenvalues.

For contrast, Hatano-Nelson has the special exact relation G_nm^OBC=e^{h(n-m)-gamma t} U_nm(t), with U unitary. Its gain depends on bounded net displacement; reverse travel cancels directional gain. There is no such established factorization for the general BdG chain.

At phi=pi/2 our chain has an exact solution. Its spatial hopping is purely imaginary, so the particle and hole hopping blocks are equal and commute with the onsite Nambu block:

$$
\omega_{m,\pm}=2J\cos\frac{m\pi}{N+1}
\pm\sqrt{\tilde\Delta^2-|g_s|^2}-i\gamma.
$$

OBC and PBC therefore have the same maximum growth at this phase. Directional drift by itself does not guarantee an OBC advantage.

## 7. Why SVD topology is not a reflection criterion

At a specified spectral reference omega_0, define

$$
Q_0=\omega_0 I-\mathcal H,\qquad
\mathcal H_{\rm ext}=\begin{pmatrix}0&Q_0\\Q_0^\dagger&0\end{pmatrix}.
$$

This extended Hermitian matrix anticommutes with diag(I,-I); its energies are plus/minus the singular values of Q_0. Small singular values and spatially separated singular vectors can give a large directional resolvent response. This does not fix the poles of the response over complex frequency, and does not bound particle-hole conversion at a boundary.

In a channel representation, a schematic round-trip operator is

$$
\mathcal L(\omega)=R_L(\omega)P_{\leftarrow}(\omega)
R_R(\omega)P_{\rightarrow}(\omega).
$$

Feedback modes satisfy det[I-L(omega)]=0; a root in Im(omega)>0 signals instability. This representation presumes a well-defined complete channel matching construction and is not evaluated in this sweep. SVD topology alone does not exclude these roots. No invariant or singular-vector topology has been calculated for the present parameter sets, and none is used to label them.

## Sources and status

- Longhi, *Convective and absolute PT-symmetry breaking in tight-binding lattices*, PRA 88, 052102 (2013), [Sec. III A and Eqs. (22)-(26)](https://arxiv.org/pdf/1310.5004#page=5): moving-ray and complex-saddle method. The paper's explicit model is Rice-Mele, not this BdG chain.
- Rassaert et al., *Emerging Non-Hermitian Topology in a Chiral Driven-Dissipative Bose-Hubbard Model*, [quadratic model, Eq. (3)](https://arxiv.org/html/2411.08965v1): model motivation and pairing convention. This document does not claim a fresh review of the published supplemental material or reproduction of Gaussian results.
- The saddle polynomial, numerical phase sweep, and finite-boundary comparison here are our calculations. A supported example is not a theorem of boundary-independent stability.
