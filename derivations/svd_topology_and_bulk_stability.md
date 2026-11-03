# SVD topology versus the bulk stability criterion

Calculated 2026-11-03 for the homogeneous quadratic background J=1, tilde Delta=0.5, |g_s|=0.4, g_c=0. This extends the phase sweep in `docs/bdg_phi_stability_results.md`; it is not a self-consistent mean-field phase diagram. The same bulk coefficients and open-chain truncation are used throughout.

## 1. Definition and source conventions

We explicitly consulted Section VI of the supplemental text included in [arXiv:2411.08965v1](https://arxiv.org/html/2411.08965v1), equations (15)-(22), in addition to the main-text winding definition (4). That version discusses the SVD of the dynamical matrix, its relation to the extended chiral Hamiltonian, and zero-frequency response. This is a review of the accessible arXiv supplemental version, not a claim to have checked a separate published supplemental file.

We retain our previously stated hopping convention: forward bond J exp(+i phi), Fourier mode exp(ikj), k integrated from -pi to pi, and time evolution exp(-i omega t). This fixes the winding sign. Reversing the spatial orientation changes its sign but not whether it is nonzero. The BdG onsite pairing is g_s, twice the Hamiltonian coefficient g_j in Rassaert Eq. (3). Amplitude loss is gamma=kappa/2. PBC closing phase is +phi, total twist N phi modulo 2pi.

At spectral reference omega_0, define

$$
Q(k;\omega_0)=\omega_0 I-\mathcal H(k),\qquad
\mathcal H_{\mathrm{ext}}(k;\omega_0)=
\begin{pmatrix}0&Q(k;\omega_0)\\Q(k;\omega_0)^\dagger&0\end{pmatrix}.
$$

The extended Hamiltonian anticommutes with diag(I,-I). Its eigenvalues are plus/minus the singular values of Q. When the bulk singular gap remains open, the integer invariant is

$$
\nu(\omega_0)=\frac{1}{2\pi i}\int_{-\pi}^{\pi}dk\,
\partial_k\log\det Q(k;\omega_0).
$$

Using -H instead of H at omega_0=0 does not change the determinant winding for this two-component problem. We focus on real omega_0=0, then check a nonzero reference. This is the homogeneous bulk invariant, so no slowly varying local-bulk approximation is invoked.

## 2. Exact zero-frequency criterion

Let

$$
a(k)=2J\sin\phi\sin k,\quad d(k)=\tilde\Delta+2J\cos\phi\cos k,
\quad \mathcal H(k)=(a-i\gamma)I+
\begin{pmatrix}d&g_s\\-g_s^*&-d\end{pmatrix}.
$$

Then

$$
\det Q(k;0)=a(k)^2-d(k)^2+|g_s|^2-\gamma^2-2i\gamma a(k).
$$

For gamma>0 and 0<phi<pi/2, the imaginary part is strictly negative for 0<k<pi and strictly positive for -pi<k<0. The curve crosses the real axis only at k=0,pi. Their real values are

$$
X_0=|g_s|^2-\gamma^2-(\tilde\Delta+2J\cos\phi)^2,
\qquad
X_\pi=|g_s|^2-\gamma^2-(\tilde\Delta-2J\cos\phi)^2.
$$

For tilde Delta>0, J>0 and the phase interval considered here, X_0<X_pi. The determinant winds once counterclockwise precisely when X_0<0<X_pi. Thus

$$
\boxed{
(\tilde\Delta-2J\cos\phi)^2+\gamma^2
<|g_s|^2<
(\tilde\Delta+2J\cos\phi)^2+\gamma^2
\quad\Longleftrightarrow\quad \nu(0)=+1.
}
$$

The invariant is zero elsewhere when the point gap is open, and undefined at gap closure. With q=sqrt(|g_s|^2-gamma^2) real, the condition can also be written

$$
|\tilde\Delta-q|<2J\cos\phi<\tilde\Delta+q.
$$

There is no nontrivial zero-frequency interval if |g_s|<=gamma. At phi=0 in our examples the determinant vanishes at interior real momenta, so the invariant is undefined rather than zero. The numerical implementation handles this endpoint explicitly instead of trusting a finite grid to detect its zeros.

## 3. Topological intervals for the tested chain

| kappa | gamma | Open phase interval with nu(0)=+1 |
|---:|---:|---|
|0.7|0.35|69.706785 degrees < phi < 81.189008 degrees|
|0.2|0.10|63.663055 degrees < phi < 86.769624 degrees|

The endpoints close the zero-frequency point gap. For other phases strictly between 0 and 90 degrees the invariant is zero; phi=90 is also gapped and trivial for these parameters. Winding sums on 8192 and 16384 momentum points agree with the exact criterion at every phase of the earlier sweep. Minimum Bloch singular values and maximum phase increments are saved as diagnostic quantities. Values very close to a transition require grid refinement or the exact criterion.

## 4. Direct comparison with stability

Rates are maximum Im(omega), in J units. The OBC values below use N=80; bulk rates are the saddle predictions, not PBC growth rates.

| phi | kappa | nu(0) | Bulk growth | OBC growth |
|---:|---:|---:|---:|---:|
|60 degrees|0.7|0|-0.045436|-0.046817|
|70 degrees|0.7|+1|-0.139335|-0.114206|
|80 degrees|0.7|+1|-0.289950|-0.168413|
|85 degrees|0.7|0|-0.350000|-0.234120|
|60 degrees|0.2|0|+0.204564|+0.203183|
|70 degrees|0.2|+1|+0.110665|+0.135794|
|80 degrees|0.2|+1|-0.039950|+0.081587|
|85 degrees|0.2|+1|-0.100000|+0.015880|

The last two rows give direct counterexamples to the proposal that zero-frequency SVD topology guarantees validity of the bulk stability classification. At phi=80, kappa=0.2, an N=40 time evolution independently confirmed growth (0.082761 fitted versus eigenvalue 0.082696). Conversely, the original phi=60 agreement occurs with nu(0)=0. At phi=85,kappa=0.7 there is substantial rate disagreement even though nu(0)=0. Distinguish agreement of the stable/unstable label from agreement of growth rates.

There is a useful exact structural observation: uniform loss shifts both OBC and bulk saddle rates by -gamma, so their difference

$$
\Delta\Gamma_N=\Gamma_{\mathrm{OBC},N}-\Gamma_{\mathrm{bulk}}
$$

is independent of gamma. In contrast, nu(0) and its transition phases depend on gamma. Therefore the zero-frequency invariant cannot universally determine whether those two rates agree. The stable/unstable label can still change as gamma crosses either rate.

In this particular sweep the nontrivial intervals overlap the region of rate disagreement. This is a correlation in a restricted slice, not an established causal explanation. No boundary channel-conversion mechanism has yet been derived, nor does the computation establish that all topology is detrimental to agreement.

## 5. Finite-OBC SVD check

Compute Q_N=-H_OBC and Q_N=U S V^dagger. The response at zero frequency is V S^-1 U^dagger when the inverse exists. This algebraic response is not a stable steady-state susceptibility if H has growing modes.

At phi=80 degrees, kappa=0.2:

| N | smallest singular value |
|---:|---:|
|20|0.01457012|
|40|0.01341337|
|80|0.00058471|
|120|0.00012521|

At N=120 the next singular value is 0.09605. The right and left singular-vector site centers are 7.97 and 113.03; about 80.7 percent of their weights lie in the first and last 12 sites respectively. Opposite-end weights are about 9.2e-6. The relation Qv=s_min u has residual 1.2e-15. This gives the expected separated singular vectors and decreasing small singular value in the same example that is dynamically unstable under OBC. Finite-size oscillations mean monotonic or asymptotic exponential scaling should not be inferred from only four sizes.

At kappa=0.7, phi=80 degrees, s_min decreases from 0.09321 (N=20) to 0.01235 (N=120); the edge localization is weaker because the parameters are close to a point-gap boundary. Nontrivial winding, rather than merely a small finite-size singular value, identifies the topology.

## 6. Reference frequency matters

The original phi=60-degree example is trivial at omega_0=0 but has nu(1.5J)=+1 for both tested losses. For kappa=0.7 at this nonzero reference, s_min is 0.02632, 0.006910 and 0.0005825 at N=20,40,80. At N=80 the right/left singular-vector centers are 7.96 and 73.04.

Thus the successful bulk prediction is not an example of a chain with no point-gap topology at any frequency. Any proposed topology/stability relation must specify the reference frequency or frequency range; a single reference cannot determine all complex-frequency poles and boundary feedback channels.

## 7. Reproduction and limits

Reusable functions: `src/bdg_topology.py`. Driver: `src/compare_topology_stability.py`. Dataset: `data/generated/bdg_topology_stability.json`. Figure: `figures/bdg_topology_stability.png`. Notebook: `notebooks/bdg_topology_stability.ipynb`. Run the driver with the same project-local environment documented in `docs/generated_data.md`.

The script checks analytic and numerical winding agreement on the complete existing phase grid, refines the momentum grid, and records finite-OBC singular values and vector localization. SVD computations use NumPy/SciPy double precision. No mean-field background is recalculated. N=80 spectra are used for comparison to avoid the strongly conditioned N=160 eigenvalues around 60 degrees; see the earlier results report for the remaining numerical limitations. The present supplement reading is explicitly the arXiv v1 Section VI; sign conventions are defined above rather than inferred from a local inhomogeneous invariant.
