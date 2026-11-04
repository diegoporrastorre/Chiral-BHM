# Literature review of boundary dependent stability in dissipative bosonic BdG chains

Literature search date: 4 November 2026. Scope: the research questions in AGENTS.md and the prescribed quadratic-chain calculations currently in derivations/. This is an extensive targeted review with reference tracing, not a guarantee that every relevant publication or priority claim has been identified.

The central phenomenon is established: a bosonic chain can be stable under OBC and unstable under PBC. Detuning-induced coupling of oppositely amplified channels, finite-size stability changes, nonnormal transient amplification, and stabilization through engineered boundary profiles also have substantial precedents. Particularly important is arXiv:2207.13728v4, which already combines nonlinear mean-field solutions, linearized fluctuations, cross-Kerr couplings, and engineered buffers. The most defensible remaining research direction is a quantitative theory for the particular density–density model and its boundary profiles, with local impulse growth, finite-chain eigenvalue stability, and nonlinear background changes distinguished explicitly.

## Scope and strength of the evidence

Searches covered bosonic Bogoliubov chains, dissipative Hatano–Nelson models, parametric resonator arrays, dynamical and topological metastability, convective versus absolute instability, non-Bloch band theory, SVD topology, chiral Bose–Hubbard models, cross-Kerr interactions, and boundary stabilization. Primary sources are linked below. Search results and aggregators were used to locate papers; substantive conclusions use author manuscripts or publisher records.

Relevant full-text sections were inspected for the closest theoretical papers. Other entries are explicitly described as abstract-level or bibliographic screening. Preprints are distinguished from journal articles. The published Rassaert main paper was checked against the local PDF. Its final publisher supplement was inaccessible through the web tool; the separately available supplementary sections in arXiv:2411.08965v1 were consulted, especially III and VI. Those sections cannot be assumed identical to the final PRL supplement.

No new numerical calculations were performed for this review. Project-specific thresholds below are existing project results, not independently revalidated here. Original reference PDFs and existing derivations were preserved.

## Direct precedents

### McDonald, Pereg-Barnea and Clerk 2018

[Phase-Dependent Chiral Transport and Effective Non-Hermitian Dynamics in a Bosonic Kitaev-Majorana Chain, Physical Review X 8, 041031](https://doi.org/10.1103/PhysRevX.8.041031); [author manuscript](https://arxiv.org/abs/1805.12557).

Sections VI–VII explicitly establish stable OBC and unstable PBC dynamics, explain circulation and repeated amplification in a ring, and analyze on-site perturbations. Boundary reflection preserves a quadrature in their ideal chain, whereas detuning converts between quadratures and can create unstable multiple-scattering paths. Their closed-chain notion of stability includes real frequencies without decay; loss must be included before comparing to strict asymptotic stability. This is a direct precedent for both the basic boundary effect and the qualitative conversion-feedback mechanism.

### Flynn, Cobanera and Viola 2020

[Deconstructing effective non-Hermitian dynamics in quadratic bosonic Hamiltonians](https://arxiv.org/abs/2003.03405), New Journal of Physics, [DOI 10.1088/1367-2630/ab9e87](https://doi.org/10.1088/1367-2630/ab9e87).

Section 4 analyzes the bosonic Kitaev chain with open, periodic, antiperiodic, and twisted boundary conditions. Krein stability and exceptional points provide an analytical treatment of boundary-sensitive spectra. This is especially relevant to closing-bond phase conventions: a ring closure is more than adding an unspecified matrix element. The work is primarily Hamiltonian, so damping conventions must be supplied in a dissipative comparison.

### Yokomizo and Murakami 2021

[Non-Bloch band theory in bosonic Bogoliubov–de Gennes systems, Physical Review B 103, 165123](https://doi.org/10.1103/PhysRevB.103.165123); [full author manuscript](https://arxiv.org/html/2012.00439v2).

Develops the generalized Brillouin zone for bosonic BdG systems and studies the Kitaev-Majorana chain, including infinitesimal-perturbation instability and reentrant behavior. It supplies a direct method for the thermodynamic OBC spectrum. A generalized-Brillouin-zone condition and a fixed-position impulse saddle condition solve different problems; their thresholds should not be equated without a model-specific argument.

### Flynn, Cobanera and Viola 2021 and 2023

[Topology by Dissipation: Majorana Bosons in Metastable Quadratic Markovian Dynamics, Physical Review Letters 127, 245701](https://doi.org/10.1103/PhysRevLett.127.245701); [manuscript](https://arxiv.org/abs/2104.03985).

Connects stable finite OBC dynamics and unstable boundaryless dynamics to long transients and topological metastability. Pseudospectra explain how asymptotically decaying systems can display substantial transient amplification. Thus a stable eigenvalue spectrum need not imply weak response or rapid relaxation.

[Topological zero modes and edge symmetries of metastable Markovian bosonic systems, Physical Review B 108, 214312](https://doi.org/10.1103/PhysRevB.108.214312); [full manuscript](https://arxiv.org/html/2306.13711v2).

Extends this framework to bosonic edge modes, approximate symmetries, correlation signatures, and nonnormal quadratic Lindbladians. Its topological metastability should be compared carefully with the project's frequency-resolved SVD amplification invariant. Shared language does not establish equality of the diagnostics or of their physical claims.

### Okuma 2022

[Boundary-dependent dynamical instability of bosonic Green's function: Dissipative Bogoliubov–de Gennes Hamiltonian and its application to non-Hermitian skin effect, Physical Review B 105, 224301](https://doi.org/10.1103/PhysRevB.105.224301); [full manuscript](https://arxiv.org/html/2202.07684v2).

Section III explicitly constructs a dissipative BdG model with instability under PBC and its absence under OBC. A dissipative Hatano–Nelson normal block competes with anomalous pairing. The paper also formulates the appropriate particle–hole symmetry through the retarded Green function. This is the closest title-level precedent for the general question asked here. Its microscopic normal block differs from a reciprocal Bose–Hubbard chain whose effective nonreciprocity emerges after linearization.

### Gómez-León, Ramos, González-Tudela and Porras 2023

[Driven-dissipative topological phases in parametric resonator arrays, Quantum 7, 1016](https://doi.org/10.22331/q-2023-05-23-1016); [full manuscript](https://arxiv.org/html/2207.13715v3).

This is the prescribed quadratic model already used in the project. Section 3.2.1 and Fig. 6 study OBC stability versus detuning, loss, and size, and explain the OBC/PBC difference using finite propagation versus recirculation. Sections 3.1 and Appendix B construct the doubled-Hamiltonian invariant; Appendices C–D use Green-function decimation. Therefore neither the model's stable topological region nor a detuning/loss/size stability study is new in itself. A high-precision recurrence or a careful separation of local impulse and finite-chain thresholds would be a more specific contribution. Published finite-size curves should be compared with the exact decoupled limit before interpreting differences as new physics or numerical errors.

### Ughrelidze, Flynn, Cobanera and Viola 2024

[Interplay of finite- and infinite-size stability in quadratic bosonic Lindbladians, Physical Review A 110, 032207](https://doi.org/10.1103/PhysRevA.110.032207); [full manuscript](https://arxiv.org/html/2405.08873v1).

The principal general-framework reference. It distinguishes finite OBC, PBC, semi-infinite, and bi-infinite systems, classifies finite/infinite stability disagreements, and treats nonnormality using pseudospectra. Coupled Hatano–Nelson and Kitaev-coupled oscillator models exhibit different types of dynamical metastability. Appendix D discusses bulk-stable but boundary-unstable systems. Its operator stability definitions are not simply fixed-site impulse decay. The project's local-decay/OBC-growth examples cannot be labeled type II solely from a saddle calculation.

## Nonlinear backgrounds and engineered boundaries

### Ramos and collaborators 2022–2024

[Topological Josephson parametric amplifier array: A proposal for directional, broadband, and low-noise amplification, arXiv:2207.13728v4](https://arxiv.org/html/2207.13728v4), revised 14 November 2024.

An essential close precedent. Equation (9) includes density–density cross-Kerr interactions together with pair hopping and number-dependent hopping. Section III derives nonlinear mean-field displacements and then linearized fluctuation dynamics. Section III.2 engineers boundary buffers through pump and detuning profiles to stabilize quasi-homogeneous amplification. Appendix D explains the mismatch between homogeneous OBC and bulk mean-field equations. This already investigates the broad methodology and stabilization strategy proposed here. The project differs through independently chosen pure density–density V, the explicit U n² convention, and its intended systematic phase diagram and boundary-scaling analysis. These differences require quantitative demonstration before supporting a novelty claim.

### Related public patent application 2025

[A method for stabilizing nonlinear parametric amplifier arrays and a superconducting quantum parametric amplifier array, WO2025196051A1](https://patents.google.com/patent/WO2025196051A1/en), published 25 September 2025; related [EP4622104A1](https://patents.google.com/patent/EP4622104A1/en).

The disclosed technical method includes Kerr and cross-Kerr resonators, nonlinear mean-field equations, linearized fluctuations, and inhomogeneous boundary controls. Figures 3, 5, and 6 and the stabilization procedure are directly relevant. It is a public technical disclosure, not a peer-reviewed universal stability theorem. Statements there about the necessity of inhomogeneity must remain tied to the modeled device and operating regime. No legal interpretation of patent scope is made here.

### Rassaert, Ramos, Roscilde and Porras 2025

[Emerging Non-Hermitian Topology in a Chiral-Driven-Dissipative Bose-Hubbard Model, Physical Review Letters 135, 203603](https://doi.org/10.1103/w32k-x3l5); [available preprint including supplementary sections](https://arxiv.org/html/2411.08965v1).

The final main paper already establishes stabilization by smoothing coherent drive and detuning near boundaries. Its steady-state calculation uses a Gaussian ansatz. Available preprint supplement III studies time evolution from vacuum and acknowledges initial-condition dependence; VI develops zero-frequency SVD and Hermitian doubling. This directly precedes profile stabilization in the target model, but does not establish the requested pure mean-field plus V phase diagram. The preprint's mean-field benchmark is not a substitute for that full diagram. The final supplement remains a version-specific verification gap.

### Dahan, Arwas and Grosfeld 2022

[Classical and quantum chaos in chirally-driven, dissipative Bose-Hubbard systems, npj Quantum Information 8, 14](https://doi.org/10.1038/s41534-022-00518-2).

Studies a small chirally driven ring using mean-field equations, fixed-point Bogoliubov stability, and quantum dynamics. It identifies stable branches, bistability, limit cycles, and chaotic attractors. This is an important PBC nonlinear precedent. A trajectory failing to settle is not by itself a demonstration that every stationary branch is unstable; existence and branch stability require separate analysis.

### Coppo, Le Boité, Felicetti and Brosco 2025–2026

[Nonlinearity-driven topology via spontaneous symmetry breaking, Physical Review B 113, 155415 (2026)](https://doi.org/10.1103/q2br-g63y); [2025 preprint](https://arxiv.org/html/2503.11928v1).

Cross-Kerr-coupled parametrically driven resonators develop topology after symmetry breaking. The preprint explicitly derives different OBC and PBC semiclassical backgrounds and resulting quadratic models. This is a close conceptual precedent for separating background reconstruction from spectral boundary effects. The model has no bare quadratic tunneling and uses a different drive and topology, so it is not the same chiral lossy Bose–Hubbard problem. Final journal and preprint details should be matched before borrowing formulas.

### Flynn, Viola and Flebus 2026

[Dynamical metastability and transient topological magnons in interacting driven-dissipative magnetic systems, arXiv:2602.13390](https://arxiv.org/html/2602.13390v1).

Extends dynamical metastability into nonlinear interacting magnetic dynamics, connecting quadratic magnon predictions to semiclassical spin evolution and transient edge physics. It prevents a broad claim that extending boundary-sensitive metastability beyond quadratic systems is untouched. Its spin platform and nonlinear constraints differ from Kerr resonators.

## Convective growth, perturbations, and long transients

### Longhi 2013 and 2016

[Convective and absolute PT symmetry breaking in tight-binding lattices, Physical Review A 88, 052102 (2013)](https://doi.org/10.1103/PhysRevA.88.052102); [manuscript](https://arxiv.org/abs/1310.5004).

Distinguishes growth seen at a fixed position from growth in a moving frame and shows how finite gain/loss regions can alter the outcome. This is the conceptual precedent for the project's impulse-response analysis. The local file named Longhi2021.pdf is this 2013 paper: its arXiv identifier is 1310.5004, despite a later internal date. Bibliographic citations should use the actual publication year.

[Quantum decay and amplification in a non-Hermitian unstable continuum (2016)](https://arxiv.org/abs/1605.02844).

Analyzes a localized state coupled to a tight-binding continuum with asymmetric hopping and distinguishes convective and absolute instability. Useful for resolvent poles and the failure of a simple correspondence between growing Bloch modes and local response. Screened at abstract level.

### Li, Lee, Mu and Gong 2020

[Critical non-Hermitian skin effect, Nature Communications 11, 5491](https://doi.org/10.1038/s41467-020-18917-4); [accessible full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC7603343/).

Shows singular effects of weak coupling between chains with different skin localization and noncommuting weak-coupling/large-size limits. This is a strong analogy for the project's detuning coupling of opposite directed channels. Establishing an exact mapping requires matching the full matrix symbol, rather than identifying only a similar spectral plot.

### Haga and collaborators 2021 and Lee and collaborators 2023

[Liouvillian Skin Effect: Slowing Down of Relaxation Processes without Gap Closing, Physical Review Letters 127, 070402](https://doi.org/10.1103/PhysRevLett.127.070402); [manuscript](https://arxiv.org/abs/2005.00824).

Demonstrates that relaxation can become slow with growing size without a closing Liouvillian gap. Relevant to interpreting long integrations, although its Liouvillian setting should not be identified automatically with a nonlinear mean-field Jacobian.

[Anomalously large relaxation times in dissipative lattice models beyond the non-Hermitian skin effect, Physical Review B 108, 064311](https://doi.org/10.1103/PhysRevB.108.064311); [manuscript](https://arxiv.org/abs/2210.14212).

Shows that propagation and interference, rather than skin localization length alone, control local relaxation. This cautions against inferring dynamics from eigenvector localization. Both entries were screened through primary abstracts and bibliographic records.

### McDonald, Hanai and Clerk 2022

[Nonequilibrium stationary states of quantum non-Hermitian lattice models, Physical Review B 105, 064302](https://doi.org/10.1103/PhysRevB.105.064302); [manuscript](https://arxiv.org/abs/2103.01941).

Shows that boundary-dependent quantum stationary densities cannot generally be inferred from the non-Hermitian drift spectrum or skin effect alone. This concerns noise-supported quantum steady states and should not be used to add covariance feedback to the project's mean-field solver. Screened at abstract level.

### Longhi 2025

[Modulational Instability and Dynamical Growth Blockade in the Nonlinear Hatano–Nelson Model](https://doi.org/10.1002/apxr.202400154); [manuscript](https://arxiv.org/abs/2501.01226).

Studies nonlinear plane-wave instability and growth blockade on a periodic Hatano–Nelson lattice. Relevant nonlinear comparison, but it is not a proof of stable stationary topological BdG amplification under OBC or a solution of the target coherently driven Kerr model.

## Topological formalism and additional evidence

The following primary sources delimit related questions. Their results should not be combined into an unrestricted claim that topology guarantees stability.

| Source | Relevance and scope |
|---|---|
| Porras and Fernández-Lorenzo, [Topological Amplification in Photonic Lattices, PRL 122, 143901 (2019)](https://doi.org/10.1103/PhysRevLett.122.143901) | SVD/Hermitian doubling and stable topological amplification in a driven photonic lattice. [Primary abstract](https://arxiv.org/abs/1812.01348) checked. |
| Herviou, Bardarson and Regnault, [PRA 99, 052118 (2019)](https://doi.org/10.1103/PhysRevA.99.052118) | SVD formulation of non-Hermitian bulk-edge correspondence; not a general dynamical stability criterion. Primary abstract checked. |
| Wanjura, Brunelli and Nunnenkamp, [Nature Communications 11, 3149 (2020)](https://doi.org/10.1038/s41467-020-16863-9) | Topological framework for directional amplification, including bosonic Kitaev quadratures and dissipation. Relevant full-text passages checked. |
| Ramos, García-Ripoll and Porras, [Topological input-output theory for directional amplification, PRA 103, 033513 (2021)](https://doi.org/10.1103/PhysRevA.103.033513) | Frequency-dependent response and the SVD mapping. Primary abstract checked. |
| Brunelli, Wanjura and Nunnenkamp, [Restoration of the non-Hermitian bulk-boundary correspondence via topological amplification, SciPost Physics 15, 173 (2023)](https://doi.org/10.21468/SciPostPhys.15.4.173) | Separates singular-vector amplification physics from OBC eigenvalue stability. Its main theorem concerns a single complex band without symmetries; applying it to a full BdG matrix requires care. Relevant full-text passages checked. |
| Peano, Houde, Marquardt and Clerk, [Topological Quantum Fluctuations and Traveling Wave Amplifiers, PRX 6, 041026 (2016)](https://doi.org/10.1103/PhysRevX.6.041026) | Stable bulk and parametrically unstable chiral edge modes in a different geometry. Important distinction between edge-mode amplification and stable directional response. Primary abstract and stability discussion checked. |
| Ling and Kain, [Selection rule for topological amplifiers in Bogoliubov–de Gennes systems, PRA 104, 013305 (2021)](https://doi.org/10.1103/PhysRevA.104.013305) | Weak-pairing stability theorem away from zero energy, applied to a different model. Not an arbitrary-loss or arbitrary-boundary theorem. Primary abstract checked. |
| Bardyn, Karzig, Refael and Liew, [Chiral Bogoliubov excitations in nonlinear bosonic systems, PRB 93, 020502(R) (2016)](https://doi.org/10.1103/PhysRevB.93.020502) | Stationary nonlinear backgrounds inducing topological Bogoliubov bands in two dimensions. Establishes broader mean-field/BdG topology precedent. Primary abstract checked. |
| Slim et al., [Optomechanical realization of the bosonic Kitaev chain, Nature 627, 767–771 (2024)](https://doi.org/10.1038/s41586-024-07174-w) | Experimental chiral amplification, size-dependent gain, and strong boundary sensitivity. Primary author/publisher abstracts checked. |
| Busnaina et al., [Quantum simulation of the bosonic Kitaev chain, Nature Communications 15, 3065 (2024)](https://doi.org/10.1038/s41467-024-47186-8) | Circuit experiment comparing open and periodic configurations and their spectra. Relevant full-text passages checked. |
| Fortin, Wang and Pereg-Barnea, [Topological amplification of the bosonic Kitaev chain with nonuniform loss, PRB 112, 064208 (2025)](https://doi.org/10.1103/d3c1-zqtp) | Loss placement changes amplification and topological transitions; useful comparison for absorbing layers. Its loss-only ideal-chain conclusions do not prove stability of detuned coupled channels. [Full manuscript](https://arxiv.org/html/2412.09744v2) inspected selectively. |
| Vega et al., [Quantum 9, 1861 (2025)](https://doi.org/10.22331/q-2025-09-18-1861) | Stable parametric topological amplification with long-range dissipative couplings. Broader context beyond nearest-neighbor models. Primary abstract checked. |
| Lee, Jin and Clerk, [Symmetry as a route to generalized bosonic Kitaev chains, arXiv:2607.08638 (2026)](https://arxiv.org/html/2607.08638v1) | Distinguishes quadrature decoupling from sublattice symmetry. Section V.4 treats chemical-potential-induced OBC instability even with stable PBC; Appendix A analyzes large-size behavior. Relevant sections inspected. |
| Ughrelidze, Viola and Cobanera, [Decoupling band topology from criticality in bosonic systems, arXiv:2607.01334v2 (2026)](https://arxiv.org/abs/2607.01334) | Additional recent separation of topology, dynamical stability, and thermodynamic stability. Abstract-level screening; not needed for the central priority assessment. |

## What the different stability criteria mean

The following distinctions are our synthesis for interpreting the project, rather than a claim that one cited paper derives every formula below.

Use i dPsi/dt = H_BdG Psi and A = -i H_BdG. A mode evolves as exp(-i omega t), so its growth exponent is Im omega = Re lambda(A). Strict stability of a finite chain requires every exponent to be negative. A neutral eigenvalue, particularly in a defective matrix, requires separate analysis.

Four distinct objects should be reported:

1. **Bloch or PBC spectral growth:** sup over real k of max Im spec H(k), with finite-ring momenta and closing twist specified when relevant.
2. **Local impulse growth:** the large-time exponent of a fixed-position Green-function matrix element or block on the bi-infinite chain. This may be negative even when a moving packet grows.
3. **Finite OBC growth:** max Im spec H_N after specifying the actual termination and profile.
4. **Infinite boundary-operator stability:** a statement about a semi-infinite or bi-infinite operator on a specified function space. It is stronger and different from decay of a chosen local matrix element.

Taking N to infinity after extracting each finite-chain eigenvalue spectrum need not recover the spectrum of an infinite operator. Taking t to infinity at fixed position need not measure growth of a norm over all positions. These differences explain why the project's local impulse criterion can fail to predict finite OBC stability without contradicting the metastability literature.

For multichannel chains, solving D(omega,z)=0 and dD/dz=0 only lists candidate spatial double roots. The causal pinch selection, numerator cancellations, and zero/infinite spatial roots must also be checked. OBC continuum bands instead involve appropriate spatial-root modulus matching and boundary conditions. Isolated boundary roots require additional consideration. Our existing detuned rank-deficient model is precisely a case where a finite-saddle-only calculation is incomplete.

### Point-gap topology at a real reference frequency

Define Q(k)=omega_ref I-H(k), with real omega_ref, and the Hermitian doubled matrix [[0,Q],[Q dagger,0]]. Its positive eigenvalues are singular values of Q. The determinant winding is

    nu(omega_ref) = (1/(2 pi i)) integral dk d_k log det Q(k).

For a finite-dimensional continuous periodic Bloch matrix, strict Bloch stability implies nu=0 at every real reference frequency. A short proof: all eigenvalues of Q lie in the upper half-plane, where a single analytic matrix logarithm exists. Hence tr log Q(k) is continuous and periodic and its total derivative integrates to zero. Therefore nonzero full determinant winding around a real reference excludes strict Bloch stability under these assumptions. This argument is included here as an explicit deduction, not a novelty claim.

It does not prove stable OBC dynamics, nor that an absorber eliminates all unstable poles. Zero total winding can hide cancellation between sectors. Complex reference points, sector-resolved invariants, and local invariants have different implications. For inhomogeneous chains, a local bulk winding additionally assumes sufficiently slow variation; it cannot establish the global eigenvalue stability of the terminated matrix.

### Background changes and the pure cross-Kerr extension

Two comparisons must remain separate: recomputing the nonlinear stationary background for each boundary condition, and changing the termination of one prescribed bulk BdG operator. The project already recognizes this distinction; the Ramos and Coppo papers make it especially important in assessing novelty.

For the convention H_on = sum_j [Delta_j n_j + U n_j²] and H_cross = V sum_bonds n_j n_k, normal ordering gives the coherent mean-field interaction term

    [U + 2 U |alpha_j|² + V sum_(k neighbor j) |alpha_k|²] alpha_j.

Thus U n² = U a†²a² + U n must be retained when translating Kerr conventions. Linearization generates normal off-site terms V alpha_j alpha_k* and anomalous terms V alpha_j alpha_k, together with normal on-site neighbor-density shifts. Those terms connect the density–density model to parametric-chain physics, but do not make it identical to a circuit Hamiltonian containing additional pair and number-dependent hopping. Signs and any absorbed linear detuning shift must be matched before comparing parameters.

## Assessment of the current project

| Proposed claim or result | Assessment after this search |
|---|---|
| Stable OBC and unstable PBC bosonic BdG chains exist | Established, with explicit theory since 2018 and dissipative BdG examples by 2022. |
| Topological amplification can coexist with stable OBC dynamics | Established in several photonic and parametric-chain models. |
| Detuning mixes opposite amplified channels and can destabilize them | Established mechanism; a new model-specific threshold remains possible. |
| Stability can depend strongly on chain size and noncommuting limits | Established general issue; exact asymptotics for a specified symbol may still add value. |
| Mean field followed by linearized BdG, with cross-Kerr interactions and boundary buffers | Already investigated in arXiv:2207.13728v4. |
| Smooth drive and detuning stabilize the chiral Bose–Hubbard chain | Already central to Rassaert et al.; its primary closure is Gaussian. |
| Pure V n_j n_(j+1), with independently controlled V approximately U and the project's exact drive/loss conventions | No exact matching systematic study identified here. This is a qualified search result, not proof of absence. |
| Local impulse decay need not imply finite OBC stability | Consistent with established distinctions; the project's particular analytic demonstration may be useful. |
| Exact recurrence and threshold for the project's detuned Quantum 2023 parameters | Potential technical contribution; priority not established by this search. |
| Fixed-width versus size-proportional taper scaling, with a fixed retained bulk | Potentially useful quantitative contribution; no exact matching result identified. |
| A right-only absorber recovers the local impulse threshold universally | Neither literature nor the project's documented finite search establishes this. |

Existing project examples sharpen the remaining questions. derivations/quantum2023_detuned_stability.md reports local impulse threshold kappa=2, a large-chain OBC threshold approximately 4.46766, and PBC threshold approximately 5.99667 for gs=gc=J=1, Delta=0.1. These values are not restated as literature results. Their importance is the coexistence of distinct thresholds in a completely specified symbol. The paper already studies stability in this model; the candidate addition is the exact distinction and its derivation.

Likewise, the N=128–512 taper comparison in derivations/bdg_long_boundary_scaling.md distinguishes a fixed boundary width from profiles occupying a fixed fraction of the chain. The right-only search in derivations/bdg_right_absorbing_boundary.md preserves an exactly unchanged core. These controlled comparisons are more specific than the already established fact that smoothing can stabilize an amplifier.

## Research direction supported by this review

A defensible framing is: **Quantify how nonlinear background reconstruction and boundary-induced channel feedback separately set the stability of a chiral Bose–Hubbard chain with independently tunable nearest-neighbor density interactions.**

The strongest next steps are:

1. Match the pure density–density Hamiltonian term by term against arXiv:2207.13728v4 and identify which extra circuit nonlinearities are absent. Record the onsite Kerr and loss normalizations explicitly.
2. Benchmark the fixed quadratic operator against Quantum 2023, Okuma 2022, and the 2024 metastability framework, including an exactly decoupled channel limit.
3. Derive a boundary matching or resolvent pole condition for the coupled model. Compare poles with local impulse singularities, rather than using topology or a few real-frequency reflection coefficients as stability proxies.
4. Determine the loss or buffer width needed as the unmodified core grows. Test whether arbitrarily weak return conversion becomes important after sufficient amplification.
5. Solve and continue mean-field stationary branches at V=0 and then V approximately U for OBC and PBC. Report residuals and BdG spectra separately from integration outcomes.
6. Test the robustness of the resulting boundary between stable and unstable regimes to size, precision, profile width, closing twist, and stationary-solver tolerance.

A useful model for the boundary calculation is a matrix round-trip condition det[I-R_L(omega) P_-(omega) R_R(omega) P_+(omega)]=0 after the relevant channels and propagation factors have been derived. This is a proposed calculation, not a formula established for the project's chain. It must be evaluated at complex frequencies and include local boundary poles. Small reflection at selected real frequencies alone cannot exclude unstable solutions.

The broad physical explanation is already in the literature. The remaining opportunity is to establish a precise stability criterion and a controlled nonlinear phase diagram for the particular model, with quantitative boundary scaling and explicit distinctions between topology, amplification, and asymptotic decay.
