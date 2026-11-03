# Chiral Bose-Hubbard project

## Scope and modeling approach

Study the driven-dissipative chiral Bose-Hubbard chain discussed in
`literature/papers/RassaertPRL2025.pdf` (Rassaert et al., Physical Review
Letters 135, 203603, 2025), extended with nearest-neighbour cross-Kerr
interactions.

Use a **mean-field approach followed by Bogoliubov analysis**, not the
self-consistent Gaussian ansatz used in the paper:

1. Derive the mean-field equations for the coherent amplitudes
   `alpha_j = <a_j>` by factorizing higher-order moments.
2. Find their stationary solutions.
3. Linearize about each stationary solution, retaining both particle and
   hole fluctuations, and calculate the non-Hermitian Bogoliubov modes.

Do not evolve Gaussian covariance matrices as part of the stationary-state
solver or feed fluctuation corrections back into the mean-field solution.

The additional interaction is `H_cross = V sum_<j,k> n_j n_k`, with each
nearest-neighbour bond counted once. For an open chain the bonds are
`(j, j+1)` for `j = 1, ..., N-1`; a periodic chain also has the closing bond.
Keep the paper's onsite interaction convention `U n_j^2` explicit when
deriving equations; do not silently substitute a different Kerr convention.

## Main research goals

### 1. Numerical steady-state phase diagram

- First establish the mean-field phase diagram at `V = 0`, then include
  nearest-neighbour cross-Kerr coupling with `V approximately U`.
- Study both open boundary conditions (OBC) and periodic boundary
  conditions (PBC).
- For OBC, compare homogeneous chains with smoothly damped boundary
  profiles of the coherent drive amplitude and detuning. Study how the
  choice of profile affects the phase diagram and dynamical stability.
- Identify stationary branches and assess their Bogoliubov stability.
  Distinguish failure of a numerical solver from evidence that a stable
  stationary state is absent.

### 2. Topology of Bogoliubov excitations

- Identify topological regions using the paper's singular-value
  decomposition (SVD) and extended chiral Hamiltonian formalism.
- State the spectral reference point or frequency and all conventions
  used to construct the extended Hamiltonian and topological diagnostics.
- Distinguish topological characterization from dynamical stability.
  If using a local bulk invariant for an inhomogeneous chain, state and
  assess the assumption of slowly varying parameters.
- Consult the paper and its supplemental material for the precise
  formalism; do not treat supplemental material as already read merely
  because it is cited in the main PDF.

### 3. Explain differences between OBC and PBC

- Investigate why stable topological regimes can occur with OBC while
  being absent with PBC, using both numerical evidence and an analytical
  argument.
- Investigate the connection to quadratic nonreciprocal chains that have
  stable regimes under OBC but not under PBC.
- Separate changes in the self-consistent mean-field background from
  changes in the spectrum and stability caused by imposing different
  boundary conditions on the same bulk Bogoliubov operator.
- Treat the proposed explanation as a research question to establish,
  rather than a conclusion to assume.

## Folder structure

```text
AGENTS.md           Project goals and instructions for coding agents.
derivations/        Analytical calculations, mean-field and Bogoliubov
                    derivations, and arguments about topology and stability.
docs/               Project documentation, conventions, numerical methods,
                    and descriptions of workflows and results.
figures/            Generated research plots and figures.
literature/
  papers/           Reference PDFs; preserve the original documents.
  notes/            Reading notes and summaries linked to their sources.
notebooks/          Exploratory calculations, parameter studies, and
                    notebooks presenting numerical results.
src/                Reusable model, solver, stability, topology, and
                    analysis code.
tmp/                Disposable intermediate files, including PDF renders.
```

Keep reusable numerical implementations in `src/` and use notebooks to
orchestrate and explain calculations. Store analytical work in
`derivations/` and literature-specific commentary in `literature/notes/`.
Temporary files are not research deliverables. If persistent numerical
datasets become necessary, establish and document a dedicated location;
the existing `.gitignore` already excludes `data/generated/`.

## Research and implementation conventions

- Document units, parameter signs, phase and gauge conventions, site
  indexing, boundary profiles, and the Bogoliubov time-evolution convention.
  Derive the stability criterion from that convention explicitly.
- For PBC, specify the closing-bond phase and total phase twist so that
  gauge transformations and OBC/PBC comparisons remain consistent.
- Record parameters, initial conditions or continuation procedures,
  convergence tolerances, and stability thresholds with numerical results.
- Check stationary residuals and sensitivity to system size and numerical
  tolerances where needed to support phase boundaries or analytical claims.
- Distinguish reproduced results, new numerical observations, and
  analytical conclusions. Do not claim reproduction of Gaussian-ansatz
  results when only the mean-field approximation has been implemented.
