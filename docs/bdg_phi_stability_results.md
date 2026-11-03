# Phase sweep: bulk saddle prediction versus open-chain stability

Calculation date: 2026-11-03. See [step-by-step derivation](../derivations/bdg_saddle_stability_step_by_step.md). This is a prescribed quadratic background, not a stationary mean-field phase diagram or Gaussian-ansatz reproduction.

Subsequent topology calculation: [SVD topology versus bulk stability](../derivations/svd_topology_and_bulk_stability.md) computes the invariant and finite-chain singular vectors. The phi=80-degree, kappa=0.2 counterexample below is nontrivial at zero frequency. References below to topology not yet being calculated describe the original phase-sweep stage.

## Main finding

The bulk saddle criterion correctly predicts the approach of the OBC stability transition near phi=52 degrees at kappa=0.7. It is not a universal OBC criterion: at phi=80 degrees and kappa=0.2 the infinite-chain fixed-site response decays, but finite open chains are unstable. Direct time evolution confirms the latter. Thus the assumption that boundary feedback is harmless fails elsewhere on this same phase sweep. No SVD topological invariant has been computed, so no topological interpretation is asserted.

## Parameters and procedure

J=1, tilde Delta=0.5, real onsite BdG pairing g_s=0.4, g_c=0. Phase from 0 to 90 degrees in steps of 5 degrees, supplemented by 52,54,56,58,62,64 degrees. Same fixed coefficients for OBC and PBC; OBC is hard truncation. Forward hopping phase is +phi; the periodic closing bond carries the same phase, total twist N phi modulo 2pi. No boundary profile, disorder, continuation of a mean field, or covariance feedback. Rates use i dPsi/dt=H Psi, so positive Im(omega) is unstable. Stability tolerance is 1e-8.

Compute the zero-loss threshold once and subtract gamma=kappa/2 exactly. Compare N=40,80,120,160 with the complex-saddle prediction. PBC curves use continuous momentum; finite rings have momentum discretization corrections.

The algebraic saddle derivative residual tolerance is 1e-6; frequency duplicates are removed within 1e-6. Positive-imaginary candidates (above 1e-8) are tested by spatial-root continuation from Omega_s+i epsilon to large positive imaginary frequency. Two runs use 1000 steps with epsilon=1e-7 and 2400 steps with epsilon=1e-8. All selected thresholds agree between these refinements. At phi=0 and 90 degrees, use exact formulas to avoid polynomial degeneracies. This is a numerical pinch diagnostic with the path limitations documented in the derivation.

## Comparison of critical amplitude-loss rates

The entries below are lossless growth rates: each is the amplitude loss gamma needed to cancel that particular rate. Strict stability requires a larger gamma. Double the entries for a critical kappa.

| phi (degrees) | Bulk saddle | OBC N=40 | OBC N=80 | OBC N=160 |
|---:|---:|---:|---:|---:|
|0|0.400000|0.398876|0.399998|0.399639|
|20|0.397993|0.378592|0.392532|0.396536|
|40|0.384184|0.372095|0.380844|0.383303|
|50|0.359144|0.351274|0.356991|0.358580|
|52|0.351009|0.343759|0.349065|0.350502|
|54|0.341558|0.335305|0.339812|0.341104|
|60|0.304564|0.301061|0.303183|0.304223*|
|65|0.263042|0.262809|0.266518|0.269050*|
|70|0.210665|0.230460|0.235794|0.238309|
|75|0.144906|0.205663|0.208837|0.211624|
|80|0.060050|0.182696|0.181587|0.185465|
|85|0.000000|0.107808|0.115880|0.116381|
|90|0.000000|0.000000|0.000000|0.000000|

An asterisk marks a dominant eigenvalue condition number above 1e9. This is a warning threshold, not a certification of other values. N=160 is especially sensitive around 58-65 degrees. Crucially, the disagreement at 75-85 degrees is already present at N=40, where conditioning is moderate, so it cannot be dismissed as large-N roundoff. Nonmonotonic finite-size convergence, including the phi=0 standing-wave sampling, is possible.

## Stability crossing at kappa=0.7

With gamma=0.35, the bulk prediction changes from growth to decay at

$$
\phi_{\rm bulk}=52.227912^\circ.
$$

The finite-OBC crossings are 50.362673 degrees (N=40), 51.779904 degrees (N=80), and 52.027622 degrees (N=120). These converge toward the saddle prediction. The Brent root solver uses angular tolerances 1e-9 radians for bulk and 1e-8 radians for OBC. Printed digits describe numerical root locations, not precision of an asymptotic boundary claim.

All sampled OBC chains and the bulk prediction are stable from 54 to 90 degrees at this loss. The continuous PBC growth remains positive until approximately 81.189 degrees (from d_min=sqrt(g^2-gamma^2)). Thus boundary-sensitive amplification is substantial even when boundaries do not change the stability classification at this particular loss.

## Counterexample: lower loss at phi=80 degrees

Keep J, tilde Delta and g_s fixed but set gamma=0.1 (kappa=0.2). The bulk saddle rate is

$$
\Gamma_{\rm abs}-\gamma=0.0600501254-0.1=-0.0399498746.
$$

| N | OBC maximum growth | Dominant eigenvalue condition | Time-domain fitted growth |
|---:|---:|---:|---:|
|40|+0.08269556|7.72|+0.08276113|
|80|+0.08158717|200.53|+0.07730682|
|120|+0.08198765|6474.42|+0.07824707|

Time propagation uses scipy.sparse.linalg.expm_multiply on -iH, starting with a real displacement at the central site: particle and hole entries both 1/sqrt(2), all others zero. Sample 501 times from 0 to 250 and fit log(total Nambu norm) over t=150..250. This computation does not use eigenvector expansion. The positive fits independently confirm instability; slower mode selection and interference explain why finite-time fits need not equal the largest eigenvalue exactly. Adding 1e-10 to the first normal diagonal shifts the maximum growth by less than 4.1e-13 in these examples; this one perturbation is a sensitivity check, not a worst-case bound.

The observation demonstrates failure of the no-boundary-instability assumption, but does not by itself identify a particular reflection-conversion path or prove an isolated boundary-mode mechanism. A transfer-matrix or boundary-resolvent study would distinguish those mechanisms.

## Independent bulk impulse check

Compute the on-site 2x2 infinite-chain propagator for delta excitations at site zero using the nonsingular matrix exponential formula. Uniform real-k trapezoidal grids have 8192 and 16384 points. Sample 601 times from 15 to 45; fit log(sqrt(t)*Frobenius norm) to account for a generic saddle prefactor.

| phi | Saddle rate | Fit t=15..45 | Range across t=15..35,25..45,15..45 |
|---:|---:|---:|---:|
|20|0.397993|0.403331|0.403331 to 0.412892|
|40|0.384184|0.383501|0.378207 to 0.383501|
|60|0.304564|0.303328|0.300934 to 0.303328|
|70|0.210665|0.210951|0.210751 to 0.211020|
|80|0.060050|0.059048|0.054867 to 0.059050|

Maximum relative quadrature change ranges from 9e-16 to 2.3e-10. Window dependence reflects oscillations and finite-time corrections; these fits support the physical saddle selection but are not high-precision estimates. In particular a finite-window slope can exceed the true asymptotic maximum. At phi=80 the propagator supports the bulk value near 0.06, not the OBC rate near 0.18.

## Artifacts and reproduction

- Figure: `figures/bdg_phi_stability.png` and vector PDF of the same name. The second panel uses kappa=0.2 and exposes the disagreement in stability classification.
- Reusable implementation: `src/bulk_stability.py`; sweep and propagation check: `src/sweep_bdg_phase.py`.
- Notebook: `notebooks/bdg_phi_sweep.ipynb`.
- Generated datasets and commands: [generated data documentation](generated_data.md).
- NumPy 2.5.4, SciPy 1.18.1, Matplotlib 3.11.2; double precision. The earlier Bloch-spectrum, uniform-loss-shift and phi=pi/2 checks remain applicable. No stationary solver is involved.

These calculations refine the earlier single-example conclusion: agreement near phi=pi/3 cannot be extrapolated to all phases, and SVD topology alone has not been established as a criterion for discarding reflections.
