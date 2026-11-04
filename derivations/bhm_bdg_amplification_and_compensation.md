# From the Bose–Hubbard model to amplification and compensated steady states

This is a pedagogical introduction to the project's mean-field and Bogoliubov approach. OBC means open boundary conditions; PBC means periodic boundary conditions (the “PBD” in the request is understood as PBC). We first explain the nonlinear background, then its linear fluctuations, then boundary-dependent stability, and finally how to engineer an exactly homogeneous background in an open chain.

The equations below are mean-field derivations. Numerical examples are explicitly identified as previously recorded project results. They do not reproduce the self-consistent Gaussian approximation of Rassaert et al.

## 1. What is the Bose–Hubbard model?

The Bose–Hubbard model describes bosons occupying a lattice of local modes. A site might represent an optical cavity or a resonator. The annihilation operator a_j removes one boson from site j, and n_j=a_j^\dagger a_j counts its occupation. Bosons hop between neighbouring sites, while interactions make the energy depend nonlinearly on occupation.

Here the chain is coherently driven and loses particles to its environment. In a frame rotating at the drive frequency, with hbar=1, take

$$
H=\sum_{j=1}^N\left[\Delta_j n_j+U n_j^2+
F_j a_j^\dagger+F_j^*a_j\right]
+J\sum_{(j,k)}\left(e^{i\phi}a_k^\dagger a_j+e^{-i\phi}a_j^\dagger a_k\right)
+V\sum_{(j,k)}n_j n_k.
$$

The detuning Delta_j is the signed coefficient displayed in this rotating-frame Hamiltonian. J>0 sets the hopping scale; U and V are signed interaction coefficients. F_j is a complex coherent drive, whose phase is experimentally meaningful relative to other drives. The baseline gauge uses real, homogeneous F_j=epsilon. All coefficients have units of frequency; one may measure them in units of J and time in units of 1/J.

Sites are labelled 1,...,N, with N>=3. Directed bonds are (j,j+1), and each interaction bond is counted once. PBC adds (N,1), with hopping J exp(i phi) from N to 1. The total ring twist is N phi modulo 2 pi. OBC has no closing bond. A gauge transformation can move hopping phases into drive phases; it must transform the drives and the ring twist consistently. Complex hopping is Hermitian here: the bare hopping magnitudes in opposite directions are equal. Directional fluctuation amplification emerges from the driven background and particle–hole dynamics, rather than from unequal bare hopping magnitudes.

Local single-particle loss is described by

$$
\dot\rho=-i[H,\rho]+\kappa\sum_j\mathcal D[a_j]\rho,
\qquad \mathcal D[a]\rho=a\rho a^\dagger-\tfrac12\{a^\dagger a,\rho\}.
$$

Thus occupation decays at rate kappa and coherent amplitude at rate gamma=kappa/2. Drive and loss mean that neither particle number nor rotating-frame energy is conserved.

## 2. The mean field is the nonlinear background

Define alpha_j=<a_j>. Its exact equation couples to higher moments. The mean-field approximation closes that hierarchy by factorizing normally ordered moments, for example <a_j^\dagger a_j a_j> approximately |alpha_j|^2 alpha_j.

Our onsite convention matters:

$$
U n_j^2=U a_j^{\dagger 2}a_j^2+U n_j.
$$

Consequently the onsite contribution to the amplitude equation is U alpha_j+2U|alpha_j|^2 alpha_j. Dropping the first term would change the stated detuning convention.

Let C be the nearest-neighbour adjacency matrix, and let T be the hopping matrix with T_{j+1,j}=J exp(i phi) and T_{j,j+1}=J exp(-i phi). Then

$$
i\dot\alpha_j=f_j=
F_j+\sum_kT_{jk}\alpha_k+
\left[\Delta_j+U-i\gamma+2U|\alpha_j|^2
+V\sum_kC_{jk}|\alpha_k|^2\right]\alpha_j.
$$

A stationary state in this frame satisfies f_j=0 at every site. This is a nonlinear complex algebraic problem: several branches may exist. A small stationary residual establishes that a candidate solves these equations; it does not establish stability or accessibility from vacuum. Likewise, a failed root search does not establish absence of a stationary branch.

In a uniform ring, alpha_j=alpha_0 and n_0=|alpha_0|^2 give

$$
0=\epsilon+\left[\Delta+U+2J\cos\phi
+2(U+V)n_0-i\gamma\right]\alpha_0,
$$
$$
\epsilon^2=n_0\left\{[\Delta+U+2J\cos\phi+2(U+V)n_0]^2+\gamma^2\right\}.
$$

The latter is a cubic in n_0. Each nonnegative real root defines a candidate uniform branch, with its phase fixed by the first equation. It does not enumerate nonuniform ring solutions.

## 3. Why the fluctuations form a BdG chain

Write a_j=alpha_j+b_j, or perturb the classical field by delta alpha_j. Expanding around a stationary background and keeping terms linear in fluctuations gives the same deterministic Bogoliubov–de Gennes (BdG) matrix. Higher-order fluctuation interactions are neglected.

The essential feature is that delta alpha and delta alpha^* mix. For example,

$$
\delta(|\alpha|^2\alpha)=2|\alpha|^2\delta\alpha+\alpha^2\delta\alpha^*.
$$

The second term converts a conjugate fluctuation into an ordinary one. In the quadratic Hamiltonian it appears as pair creation and annihilation: the driven coherent background supplies energy for parametric amplification.

At a stationary background,

$$
i\partial_t\delta\alpha=A\delta\alpha+B\delta\alpha^*,
$$
$$
A_{jk}=T_{jk}+\delta_{jk}\left[\Delta_j+U-i\gamma
+4U|\alpha_j|^2+V\sum_lC_{jl}|\alpha_l|^2\right]
+VC_{jk}\alpha_j\alpha_k^*,
$$
$$
B_{jk}=2U\alpha_j^2\delta_{jk}+VC_{jk}\alpha_j\alpha_k.
$$

Onsite interactions create onsite pairing; cross-Kerr interactions additionally modify normal intersite couplings and create intersite anomalous pairing. The doubled, or Nambu, vector and its generator are

$$
z=\begin{pmatrix}\delta\alpha\\\delta\alpha^*\end{pmatrix},
\qquad i\dot z=\mathcal Lz,
\qquad \mathcal L=\begin{pmatrix}A&B\\-B^*&-A^*\end{pmatrix}.
$$

The conjugate component is often called the hole component. It is necessary to represent amplitude and phase fluctuations; it is not an additional independent classical field. Operator Langevin equations also contain reservoir noise. Noise seeds fluctuations, but the deterministic matrix above decides their linear stability. We do not evolve covariances or feed fluctuation corrections into alpha.

Although the underlying coherent Hamiltonian is Hermitian, the bosonic dynamical matrix is generally non-Hermitian because of its particle–hole structure and loss. Its eigenvalues are frequencies, with pairing omega <-> -omega^*. For a mode z(t)=z(0) exp(-i omega t),

$$
|e^{-i\omega t}|=e^{\operatorname{Im}\omega\,t}.
$$

Therefore strict finite-chain linear asymptotic stability means max Im(omega)<0. Positive imaginary part means exponential instability. A zero value needs separate treatment of algebraic growth and Jordan blocks. Negative eigenvalue growth rates still allow substantial transient amplification when eigenvectors are nonorthogonal.

## 4. Convective amplification versus growth at a fixed position

Imagine injecting a short, localized pulse. Two different questions can be asked: does the pulse grow while travelling, and does the response continue growing at the injection site? They need not have the same answer.

For an infinite uniform BdG chain, define the impulse-response matrix

$$
G_r(t)=\int_{-\pi}^{\pi}\frac{dk}{2\pi}\,
e^{ikr}e^{-i\mathcal L(k)t},
\qquad
\lambda(v)=\limsup_{t\to\infty}\frac1t\log\|G_{\lfloor vt\rfloor}(t)\|.
$$

Here lattice spacing is one, r is displacement from the source, and v is the observation velocity. Exponential local growth is measured by lambda(0).

* **Convective instability:** lambda(0)<0 but lambda(v)>0 along some moving ray. A travelling pulse grows while each fixed observation site eventually sees it decay.
* **Absolute instability:** lambda(0)>0. A localized disturbance grows even at a fixed position.
* A stable driven amplifier can also amplify a continuously injected signal spatially. That steady forced response is a different measurement from the long-time impulse exponent.

Thus “local amplification” needs a specified meaning: finite-time gain at a site, amplification of a local drive, and asymptotic absolute instability are different statements. The convective/absolute distinction used here follows [Longhi's tight-binding analysis](https://arxiv.org/abs/1310.5004); applying it to this BdG chain requires its own dispersion and boundary analysis.

A useful continuum illustration is

$$
\partial_tu+v_d\partial_xu=\mu u+D\partial_x^2u,\qquad D>0.
$$

Its localized Green function is proportional to

$$
t^{-1/2}\exp\left[\mu t-\frac{(x-v_dt)^2}{4Dt}\right].
$$

At the pulse centre, x=v_dt, the growth rate is mu. At fixed x it is mu-v_d^2/(4D). Hence 0<mu<v_d^2/(4D) gives convective growth with local decay. This is an illustrative transport equation, not a derived continuum limit of our interacting chain.

For lattice bands omega_s(k), the analogous moving-ray calculation uses complex saddles satisfying omega_s'(k_s)=v, with rate Im[omega_s(k_s)-v k_s]. Only saddles contributing to the physical Fourier contour count; maximizing over every algebraic saddle can give a wrong answer. See [the step-by-step saddle derivation](bdg_saddle_stability_step_by_step.md).

## 5. How OBC can be stable while PBC is unstable

First hold the bulk BdG coefficients fixed. This isolates the spectral effect of changing boundaries from changes in the nonlinear stationary background.

For V=0 and a uniform background, rotate its pairing to a real g=|2U alpha_0^2| and define delta=Delta+U+4Un_0. With Fourier convention exp(ikj),

$$
\mathcal L(k)=[2J\sin\phi\sin k-i\gamma]I+
\begin{pmatrix}
d(k)&g\\-g&-d(k)
\end{pmatrix},\qquad d(k)=\delta+2J\cos\phi\cos k,
$$
$$
\omega_\pm(k)=2J\sin\phi\sin k-i\gamma
\pm\sqrt{d(k)^2-g^2}.
$$

On a ring, k=2 pi m/N. Any allowed mode with Im(omega)>0 makes the ring unstable. A periodic disturbance cannot leave the system. In the infinite-chain description, however, a localized disturbance combines many momenta and can be carried away while growing. The odd-in-k term supplies drift, and the pairing permits growth.

With open ends there is no imposed periodic return path. A convectively amplified disturbance can leave the active region or be dissipated, permitting negative growth rates for all finite-chain eigenmodes. This is the physical route to OBC stability with PBC instability. Hard open boundaries also scatter waves and mix BdG channels: they are not automatically perfect absorbers. Reflections or boundary-localized modes may create unstable feedback. Local bulk decay is therefore a diagnostic, not a general proof of finite-OBC stability. A closed-loop explanation must include propagation, reflection, and particle–hole conversion, rather than assume that any reflected wave experiences the same gain.

One previously recorded **fixed quadratic operator** example uses J=1, phi=pi/3, delta=0.5, g=0.4, gamma=0.35. The continuous-PBC maximum growth is +0.05; rings divisible by three attain it at k=+/-2 pi/3. The recorded N=40 OBC growth is approximately -0.04894, and N=80 gives -0.04682. These are existing eigenvalue results, classified with tolerance 1e-8, not new computations or nonlinear stationary solutions. The size study and its large-N conditioning limitations are in [the original comparison](convective_bdg_stability.md).

Drift alone does not guarantee this effect. At phi=pi/2, the onsite-paired chain has exact OBC frequencies

$$
\omega_{m,\pm}=2J\cos\frac{m\pi}{N+1}
\pm\sqrt{\delta^2-g^2}-i\gamma,
$$

and its maximum growth equals the PBC value. This counterexample prevents interpreting chirality as a universal stabilizer.

For the nonlinear BHM, simply cutting a ring generally changes alpha_j, and hence both A and B. A clean comparison has two parts: compare boundaries on prescribed identical bulk coefficients; then solve and compare the actual stationary backgrounds. The compensation below helps make the first comparison physically realizable. For V nonzero, edge coordination also changes diagonal BdG terms, so even a homogeneous background requires care in defining an identical-operator comparison.

Large amplification is also compatible with stable eigenvalues. A monochromatic force probes (omega I-L)^{-1}, whose norm is the inverse smallest singular value of omega I-L. A nonnormal stable matrix can have a large resolvent. This observation does not by itself establish a topological invariant, nor does a topological label establish dynamical stability.

## 6. Compensation: replace the missing mean-field forces

### 6.1 Why a homogeneous drive usually gives a nonuniform open chain

Choose a uniform bulk solution alpha_0 of the cubic above. In the interior, hopping supplies 2J cos(phi) alpha_0 and cross-Kerr interactions supply 2Vn_0 alpha_0. At an open endpoint there is only one neighbour. With all drives and detunings left uniform, the endpoint equation therefore differs from the bulk equation. Uniform density alone is not sufficient: the complex amplitudes must satisfy the phase-sensitive equations as well.

Compensation means adding known static controls that restore the missing force at the chosen target. It is branch dependent: different alpha_0 require different compensation.

### 6.2 Exact complex-drive compensation

At the left endpoint the missing hopping force is J exp(+i phi) alpha_0; at the right it is J exp(-i phi) alpha_0. Each endpoint is also missing Vn_0 alpha_0 from the absent interaction bond. Keep Delta_j=Delta and choose

$$
\boxed{F_1=\epsilon+[J e^{i\phi}+Vn_0]\alpha_0,\qquad
F_N=\epsilon+[J e^{-i\phi}+Vn_0]\alpha_0,}
$$
$$
F_j=\epsilon\quad(2\le j\le N-1).
$$

The sign is positive: the additional drive replaces the term removed by cutting the bond. To check it, the uncompensated left residual evaluated at alpha_0 is

$$
f_1^{\rm uncomp}=f_{\rm bulk}-[J e^{i\phi}+Vn_0]\alpha_0.
$$

Adding the boxed drive correction makes it exactly zero; the right endpoint works identically. No extra physical neighbour or fluctuating closing bond is introduced. For V=0 these reduce to the project's endpoint additions h_L=J exp(i phi) alpha_0 and h_R=J exp(-i phi) alpha_0.

These are complex field amplitudes, not pump intensities. In general both phase and magnitude must be controlled. The target is exactly homogeneous at mean-field level, including phase; quantum fluctuations need not be homogeneous near an edge.

### 6.3 Detuning compensation for the cross-Kerr part

An equivalent way to replace the missing cross-Kerr mean-field force is

$$
\Delta_1=\Delta_N=\Delta+Vn_0,
\qquad F_1=\epsilon+J e^{i\phi}\alpha_0,
\qquad F_N=\epsilon+J e^{-i\phi}\alpha_0.
$$

Interior parameters stay unchanged. This produces the same target stationary amplitude as drive-only compensation. However, these implementations have different BdG matrices: detuning contributes to the derivative of the force, whereas a fixed drive does not. Detuning compensation restores the missing diagonal Hartree contribution at the edge, while the removed intersite normal and anomalous couplings remain absent. It must therefore be specified when reporting stability for V nonzero.

For V=0, fixed endpoint-drive compensation has zero derivative with respect to fluctuations. Its BdG matrix at alpha_0 is precisely the hard-OBC truncation of the homogeneous bulk BdG operator. A fictitious neighbour that followed the fluctuations, or active feedback that changed with alpha_j(t), would define a different problem.

### 6.4 General inverse construction

For a chosen target alpha_j^tar, any specified detuning, hopping, interaction graph, and uniform loss can be balanced by

$$
F_j^{\rm tar}=-\sum_kT_{jk}\alpha_k^{\rm tar}
-\left[\Delta_j+U-i\gamma+2U|\alpha_j^{\rm tar}|^2
+V\sum_kC_{jk}|\alpha_k^{\rm tar}|^2\right]\alpha_j^{\rm tar}.
$$

This is a direct substitution into f=0, so it establishes existence of the target for those controls. It is also applicable to nonuniform targets or boundary profiles. Smoothly damping the drive and detuning is a separate design choice; such profiles do not generally preserve a uniform state unless the inverse equation is enforced.

## 7. What compensation guarantees, and how to use it

Compensation guarantees a stationary solution of the specified mean-field equations. Stability still requires the BdG spectrum. Preparation from vacuum requires a dynamical protocol and a suitable basin of attraction. Stability under finite perturbations is a further nonlinear question.

In practice:

1. Choose J, phi, Delta, U, V, kappa and epsilon; solve the bulk cubic and select a branch n_0. Recover its complex alpha_0.
2. Apply the branch-specific fixed endpoint controls, documenting whether cross-Kerr compensation uses drives or detunings.
3. Substitute alpha_j=alpha_0 into every stationary equation and record max_j|f_j|, its units, and the numerical tolerance.
4. Construct the full OBC BdG matrix and compare its maximum growth with that of the specified ring, using its discrete momenta and the same closing-bond phase. Report size and stability tolerance; check conditioning when margins are small.
5. Test preparation and sensitivity separately. If controls are imperfect, solve for the displaced background before assigning its actual BdG stability.

The last point matters particularly in an amplifier. For a small fixed drive error delta F,

$$
\mathcal L\begin{pmatrix}\delta\alpha\\\delta\alpha^*\end{pmatrix}
=-\begin{pmatrix}\delta F\\-\delta F^*\end{pmatrix}
$$

is the linear stationary-response equation, provided the stationary Jacobian is invertible. A small singular value at frequency zero can turn a tiny compensation error into a large background deformation even while all OBC eigenmodes decay. The nominal endpoint correction need not grow with N, while the required precision can become more demanding with length. Direction and input quadrature matter; this is not a universal exponential law for every error. See [the project's compensation sensitivity derivation](compensation_sensitivity_scaling.md) for a quantified example and its precision checks.

## Further reading and provenance

* [Mean-field and Bogoliubov equations](mean_field_bogoliubov.md): the project's full conventions and cross-Kerr Bloch matrix.
* [Bulk saddle analysis and boundary caveats](bdg_saddle_stability_step_by_step.md): mathematical details behind the impulse-response discussion.
* [Fixed-background numerical comparison](convective_bdg_stability.md): provenance and limitations of the quoted OBC/PBC example.
* [Longhi, Physical Review A 88, 052102 (2013)](https://arxiv.org/abs/1310.5004): convective versus absolute instability in tight-binding systems; its explicit model differs from this BdG chain.
* [Rassaert et al., project reference PDF](../literature/papers/RassaertPRL2025.pdf): physical motivation. This introduction uses the project's mean-field derivation and does not claim a fresh reading of the supplemental material or reproduction of its Gaussian results.
