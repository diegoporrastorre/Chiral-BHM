# Longer chains with Rassaert-inspired boundary profiles

This extends `bdg_rassaert_damped_boundaries.md` to N=128,256,512.
The prescribed quadratic model is unchanged: bulk gs=gc=J=1, Delta=0.1,
phi=pi/2, P=0, hbar=1, lattice spacing one. Uniform Lindblad loss kappa
gives amplitude damping kappa/2. We use i partial_t Psi=H Psi, so a
positive imaginary eigenfrequency means growth.

For OBC, with sites j=1,...,N, the joint taper remains

$$
s_j=\tanh\left(\frac{\min(j,N+1-j)}{N_0}\right),\qquad
\Delta_j=0.1s_j^2,\quad g_{s,j}=s_j^2,\quad
g_{c,j,j+1}=s_js_{j+1}.
$$

The hopping and loss stay uniform. This is the previously specified
adaptation of the paper's smooth profile, not a nonlinear stationary
background or an added absorbing-loss layer. No covariance evolution,
mean-field recalculation, or additional boundary scattering term is used.

## Stability thresholds

Uniform loss shifts every generator eigenvalue by -kappa/2. Thus

$$
\Gamma_N(\kappa)=\frac{\kappa_c(N)-\kappa}{2},\qquad
\text{strict stability requires }\kappa>\kappa_c(N).
$$

The following are numerical critical losses in units J. Values are
rounded to four decimals; the saved brackets are narrower than 2e-5.

| N | Fixed N0=8 | N0=N/8 | N0=N/5 |
|---:|---:|---:|---:|
| 128 | 3.9204 | 3.4289 | 2.7535 |
| 256 | 4.1932 | 3.4635 | 2.7856 |
| 512 | 4.3300 | 3.4811 | 2.8018 |

At kappa=3.5, all tested fixed-width N0=8 chains are unstable, whereas
all tested N0=N/8 and N0=N/5 chains are stable. At kappa=3, only the
N0=N/5 series is stable among these three choices. Equality at a
threshold does not give strict exponential stability.

## What changes when the boundary width scales?

With fixed N0=8, the near-uniform interior grows while the transition
regions remain of fixed width. Its increasing critical loss is consistent
with weakening protection from a fixed boundary taper; this data alone
does not prove its infinite-N limit. For comparison, the sharp-boundary
limit established earlier is 2 sqrt(4.99)=4.46766158.

For N0=N/8, the central s_j^2 is tanh(4)^2=0.99865905 at all even sizes.
The central Delta and pairings remain within 0.1341% of their original
bulk values. However, the modified region occupies a fixed fraction of
the chain: this is not a fixed-size boundary attached to an ever-longer
uniform core. Approximately one quarter of the sites have s_j^2>0.99.

For N0=N/5, the central s_j^2 is tanh(2.5)^2=0.97340777. The center is
already 2.6592% below the original bulk coefficients, and no site has
s_j^2>0.99. Its stronger stabilization must be interpreted together
with that finite change throughout the chain.

The original uniform local-bulk criterion remains kappa>2 as a reference.
None of these longer-chain thresholds restores it. The data demonstrate
stabilization through N=512 for the specified scaled profiles, not
stability at every length or an exact infinite-size threshold.

## Numerical method and verification

Ordinary long-chain double-precision diagonalization is not used to
determine these thresholds. The real site-interleaved generator has
bandwidth two. We construct its characteristic polynomial and apply
Routh counts to count ALL roots to the right of trial growth rates.
Profile coefficients themselves are evaluated with Decimal arithmetic.

Working precisions are N+80 decimal digits: 208,336,592 for these sizes.
Both endpoints of each final growth-rate bracket are checked again with
100 additional digits: 308,436,692. The higher-precision counts are positive
to the right of the lower endpoint and zero to the right of the upper
endpoint for every reported case. These checks establish numerical
dominance of the root, rather than just finding some characteristic root.
They are precision-refined numerical checks, not interval certificates.

To reduce runtime, some cases first locate a real-root candidate by
evaluating the characteristic polynomial, then verify it using full
Routh counts. The code falls back to Routh bisection if the candidate
is not the spectral abscissa. Other cached cases used Routh bisection
throughout; both methods use the same higher-precision final checks.
The maximum requested growth bracket width is 1e-5. Use the saved upper
critical-loss bracket for a conservative numerical stability decision.

No new long-chain time-domain fit is used as independent evidence here:
the independent time-evolution checks in the preceding study were at
N=40. Extending that check naively in double precision would encounter
the same nonnormal sensitivity. The longer-chain claims here rest on
the refined characteristic-polynomial root counts and exact uniform-loss
shift, not on a raw eigensolver or a short transient.

## Reproduction and files

Run `python src/bdg_long_boundaries.py` from the project root. It resumes
missing cases from `data/generated/bdg_long_boundaries.json`; remove or
rename that dataset explicitly if a full recomputation is desired.
The implementation uses `src/banded_stability_decimal.py`. All parameters
are deterministic; there are no random seeds or stationary initial guesses.

The figure is `figures/bdg_long_boundaries.png`. The notebook
`notebooks/bdg_long_boundaries.ipynb` loads the results, displays the
brackets and profiles, and provides a recomputation entry point.
These are new fixed-BdG numerical observations, not a reproduction of
the nonlinear phase diagram of Rassaert et al.
