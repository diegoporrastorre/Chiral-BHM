# Mean field and Bogoliubov equations

Set hbar = 1. All energies and decay rates can be expressed in units of J,
and time in units of 1/J. Site labels below run from 1 to N.

## Hamiltonian and boundary conventions

In the gauge with real coherent drive amplitudes,

$$
H=\sum_j[\Delta_j n_j+U n_j^2+\epsilon_j(a_j+a_j^\dagger)]
+\sum_{(j,k)}J(e^{i\phi}a_k^\dagger a_j+e^{-i\phi}a_j^\dagger a_k)
+V\sum_{(j,k)}n_jn_k.
$$

Directed bonds are (j,j+1), plus (N,1) for PBC. Each undirected bond is
counted once; we assume N >= 3. Loss is kappa times the standard Lindblad
dissipator D[a_j], so amplitudes decay at kappa/2.

Let C be the real symmetric nearest-neighbour adjacency matrix and T the
Hermitian hopping matrix: T[k,j] = J exp(i phi) on each directed bond.
The PBC model uses this phase also on the closing bond, with total twist
N phi modulo 2 pi. It is not automatically gauge equivalent to a ring
with real hopping and a periodic drive phase: that equivalence requires
a compatible twist (N phi = 0 modulo 2 pi for an untwisted ring).

For homogeneous profiles, epsilon_j = epsilon and Delta_j = Delta.
For damped OBC, let d_j = min(j,N+1-j) and s_j = tanh(d_j/N0); use
epsilon_j = epsilon s_j and Delta_j = Delta s_j^2. N0 > 0 is the boundary
width in lattice sites. This mirrors the profile in Eq. (4) of the paper.

## Mean-field equation

Normal ordering matters: n_j^2 = a_j^{dagger 2} a_j^2 + n_j.
Factorizing normally ordered moments therefore gives

$$
i\dot\alpha_j=f_j(\alpha,\alpha^*)=
\epsilon_j+\sum_kT_{jk}\alpha_k+
\left[\Delta_j+U-\frac{i\kappa}{2}
+2U|\alpha_j|^2+V\sum_k C_{jk}|\alpha_k|^2\right]\alpha_j.
$$

In the OBC bulk, the hopping contribution is
J(exp(-i phi) alpha_{j+1} + exp(i phi) alpha_{j-1}).
The +U shift must be retained for the U n^2 convention. Solve f = 0 for
stationary states. No Gaussian covariance enters these equations.

## Bogoliubov matrix

At a stationary solution alpha, write i d(delta alpha)/dt =
A delta alpha + B delta alpha*, with Wirtinger derivatives
A[j,k] = partial f_j / partial alpha_k and
B[j,k] = partial f_j / partial alpha_k*:

$$
A_{jk}=T_{jk}+\delta_{jk}\left[\Delta_j+U-i\kappa/2
+4U|\alpha_j|^2+V\sum_l C_{jl}|\alpha_l|^2\right]
+VC_{jk}\alpha_j\alpha_k^*,
\qquad
B_{jk}=2U\alpha_j^2\delta_{jk}+VC_{jk}\alpha_j\alpha_k.
$$

For the Nambu vector z = (delta alpha, delta alpha*)^T,

$$
i\dot z=\mathcal L z,\qquad
\mathcal L=\begin{pmatrix}A&B\\-B^*&-A^*\end{pmatrix}.
$$

Modes evolve as exp(-i omega t). Strict linear asymptotic stability
requires max Im(omega) < 0; a positive value is instability. Values near
zero within a declared numerical tolerance are marginal/unresolved.
The spectrum has particle-hole pairing omega <-> -omega*.

## Homogeneous PBC benchmark

A uniform solution has density n and obeys

$$
\epsilon^2=n\{[\Delta+U+2J\cos\phi+2(U+V)n]^2+(\kappa/2)^2\},
\qquad
\alpha=-\frac{\epsilon}{\Delta+U+2J\cos\phi+2(U+V)n-i\kappa/2}.
$$

Every nonnegative real root is a candidate uniform branch. This cubic
does not exclude additional nonuniform PBC solutions.
With delta alpha_j proportional to exp(i q j), q = 2 pi m/N,

$$
A(q)=\Delta+U+4Un+2Vn+2J\cos(q-\phi)+2Vn\cos q-i\kappa/2,
\quad B(q)=2\alpha^2(U+V\cos q),
$$
$$
\mathcal L(q)=\begin{pmatrix}A(q)&B(q)\\-B(q)^*&-A(-q)^*\end{pmatrix}.
$$

These formulas provide an independent check of the real-space spectrum.
Topology and the SVD/extended chiral construction will be developed
separately after checking the relevant supplemental formalism.
