# Fixed-background BdG stability: first saddle/OBC comparison

Date: 2026-11-03. New analytical application and numerical observations, not a reproduction of Rassaert's Gaussian calculation or a mean-field phase diagram.

Follow-up: the wider phase sweep is documented in `docs/bdg_phi_stability_results.md`, with a step-by-step derivation in `derivations/bdg_saddle_stability_step_by_step.md`. It finds a genuine disagreement between the bulk criterion and OBC stability at phi=80 degrees, gamma=0.1. The single-example agreement below must not be extrapolated to arbitrary phase.

## Conventions

Homogeneous quadratic chain, onsite pairing only (g_c=0). hbar and lattice spacing are one. J>0; onsite real detuning delta denotes tilde Delta. A global phase rotation makes onsite BdG pairing g real and positive. The Hamiltonian pairing term is (g b_j^dagger^2 + h.c.)/2. gamma=kappa/2 is amplitude loss. Nambu ordering is (all b, all b^dagger), with i dPsi/dt=H Psi; max Im(omega)<0 means strict stability. Sites are 1,...,N. OBC has hard truncation, no boundary damping profile. PBC closes the forward hopping with J exp(i phi), total twist N phi modulo 2pi. The same coefficients are used under both boundaries; no stationary background is recomputed.

For b_j proportional to exp(ikj),

    H(k) = [a(k)-i gamma] I + [[d(k),g],[-g,-d(k)]],
    a(k)=2 J sin(phi) sin(k), d(k)=delta+2 J cos(phi) cos(k),
    omega_±(k)=a(k)-i gamma ± sqrt(d(k)^2-g^2).

PBC growth in the continuous-k limit is

    Gamma_PBC = -gamma + sqrt(max(0,g^2-d_min^2)),
    d_min=max(0,abs(delta)-2 abs(J cos(phi))).

Finite rings require their discrete momenta. No spectral topology is calculated here.

## Applying Longhi's method

Source: Stefano Longhi, PRA 88, 052102 (2013), https://doi.org/10.1103/PhysRevA.88.052102 ; https://arxiv.org/pdf/1310.5004 . Section III A, equations (22)-(26), supplies the moving-ray/saddle method. His explicit model is Rice-Mele; confinement of gain/loss in that paper is not identical to imposing hard OBC here.

For the infinite-chain localized response, evaluate the matrix Fourier integral exp(ikn) exp[-iH(k)t]. Along n=vt, candidate saddles obey omega'(k)=v and their exponents are Im(omega-vk). Only saddles accessible to the integration contour contribute. Algebraic roots alone do not perform this selection.

The implementation constructs stationary candidates using z=exp(ik), s=J sin(phi), c=J cos(phi), q=delta z+c(z^2+1):

    s^2 (z^2+1)^2 [q^2-g^2 z^2] + c^2 q^2 (z^2-1)^2 = 0.

Squaring introduces extra roots; the original branch derivative is checked to 1e-6. Degenerate square-root points are excluded. This is not a general pinch-tracking algorithm and can fail in degenerate parameter limits. It must not be used by simply maximizing all candidate imaginary parts.

For example, delta=0, g=1.2, J=1, phi=pi/3 produces a candidate Im(omega)=1.32665, exceeding the real-k maximum 1.2. It cannot govern the physical response. This illustrates why contour selection matters.

## Numerical example

J=1, phi=pi/3, delta=0.5, g=0.4. Without loss, candidate growing saddles are

    k_s = ±1.878582121893184 + 0.191318746557743 i,
    omega_s = ±1.594344908205302 + 0.304564113599121 i.

These are supported as contributing dominant saddles by a direct on-site propagator integral, rather than by a completed analytic contour deformation. Real-k maximum growth is 0.4. The resulting candidate bulk convective interval is

    0.3045641136 < gamma < 0.4,
    equivalently 0.6091282272 < kappa < 0.8.

This is a supported bulk prediction, not a proven finite-OBC boundary for every N.

Set gamma=0.35 (kappa=0.7). PBC continuous-k growth is +0.05. k=±2pi/3 attains it, so finite rings divisible by three also attain it.

| N | OBC growth without loss | OBC growth with gamma=0.35 |
|---|---:|---:|
|20|0.2880346094|-0.0619653906|
|40|0.3010606343|-0.0489393657|
|80|0.3031826689|-0.0468173311|
|160|0.3042229874|-0.0457770126|
|240|0.3056475874|-0.0443524126|

The N=240 value breaks the convergence trend. Treat it as unresolved large-N numerical sensitivity, not evidence of a shifted physical threshold. Tiny backward residuals do not guarantee accurate eigenvalues in a nonnormal problem. Agreement is persuasive through N=160, but no uniform-size stability proof is claimed. Classification tolerance: 1e-8; all tabulated negative margins are much larger.

Direct infinite-chain impulse response: initial unit excitation at site zero in each Nambu component, assembled into the on-site 2x2 Green matrix. Integrate on uniform grids of 4096, 8192, 16384 momenta. Frobenius norm sampled at 801 times from t=20 to 60. Linear fit of log(sqrt(t)*norm) gives 0.3049786866 without loss; restricting to t=40..60 gives 0.3083269647. Oscillatory interference causes fit-window variation. Quadrature refinement changes the norm by less than 2e-13 relative over this interval. Subtract gamma=0.35 for the lossy exponent. The response supports approximately 0.3046, but the finite-time fit is not an independent high-precision determination. Longer times can suffer cancellation among growing Fourier contributions.

## Exact counterexample to 'drift alone stabilizes OBC'

At phi=pi/2 the real-space hopping matrix T is purely imaginary, so -T*=T. Consequently the lossless BdG matrix separates into commuting spatial and onsite Nambu parts:

    H = I_Nambu tensor T + [[delta,g],[-g,-delta]] tensor I_sites - i gamma I.

The OBC spatial eigenvalues are real, 2J cos(m pi/(N+1)). Thus

    omega_(m,±) = 2J cos(m pi/(N+1)) ± sqrt(delta^2-g^2) - i gamma.

Both boundaries have maximum growth sqrt(max(0,g^2-delta^2))-gamma. There is no OBC/PBC stability window at this phase. At g=abs(delta), gamma=0 the onsite block can be defective: zero exponential growth does not imply bounded evolution. Positive gamma removes that marginal issue.

## Reproducibility and checks

Implementation: src/bulk_stability.py. Notebook: notebooks/convective_bdg_stability.ipynb. NumPy 2.5.4, SciPy 1.18.1; double-precision eigensolver and uniform trapezoidal Fourier quadrature. No random seeds, root-solver initial guesses, continuation, or stationary residuals apply because this is a prescribed quadratic operator.

benchmark() checks N=20 PBC eigenvalues against Bloch dispersion (maximum matching error 1.07e-14), exact uniform-loss spectral shift (1.61e-14), and the phi=pi/2 analytical growth (2.39e-15). Dominant eigenpair relative residuals are below 7e-16, with the conditioning limitation above. No persistent numerical dataset is introduced; the table is the recorded result.

Next steps before a general phase boundary: implement pinching/contour continuation, address large-N eigenvalue conditioning, then test boundary-mode exceptions and inhomogeneous stationary backgrounds. These results establish a concrete example, not a general parameter-space classification.
