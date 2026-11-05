# Drive width near an intermediate-branch merger

For the uniform background with nearest-neighbour cross-Kerr interaction,
use the onsite convention `U n_j^2`, hopping `J=1`, and

$$
c=\Delta+U+2\cos\phi,\qquad \beta=2(U+V),\qquad
\gamma=\kappa/2.
$$

The compensated OBC solution has the same algebraic bulk equation as the
uniform PBC solution,

$$
\epsilon^2=E(n)=n[(c+\beta n)^2+\gamma^2].
$$

This equality concerns the stationary background. It does not equate the
OBC and PBC fluctuation operators: the compensated OBC endpoints retain
their degree-one Hartree derivatives.

Assume `beta < 0` and `c > sqrt(3) gamma`, and define

$$
R=\sqrt{c^2-3\gamma^2},\quad
n_-= -\frac{2c-R}{3\beta},\quad
n_+= -\frac{2c+R}{3\beta},\quad
n(f)=n_-+f(n_+-n_-).
$$

Distinct intermediate states have `0 < f < 1`. The lower-density fold
`f=0` has the upper drive `epsilon_fold=sqrt(E(n_-))`. Since `E` is cubic,
its Taylor expansion at that fold is exact:

$$
E(n_-+\delta n)=E(n_-)+\beta R\,\delta n^2+\beta^2\delta n^3,
\qquad \delta n=-\frac{2R}{3\beta}f.
$$

Consequently,

$$
\epsilon_{\mathrm{fold}}^2-\epsilon(f)^2
=-\frac{4R^3}{27\beta} f^2(3-2f),
$$

and the exact positive drive difference is

$$
\boxed{\epsilon_{\mathrm{fold}}-\epsilon(f)
=\frac{-4R^3 f^2(3-2f)}
{27\beta[\epsilon_{\mathrm{fold}}+\epsilon(f)]}.}
$$

If a stable OBC component extends from the fold to a first instability at
`f=f_c(N)`, its drive width is this expression evaluated at `f_c(N)`.
For a narrow component,

$$
\Delta\epsilon(N)
=-\frac{2R^3}{9\beta\epsilon_{\mathrm{fold}}}\,f_c(N)^2
+O(f_c(N)^3).
$$

Thus a small density interval becomes an even smaller drive interval.
The narrow drive width alone does not establish marginal dynamical
stability or a topological mechanism. The growth rate must be evaluated
separately. If the leading growth has a locally smooth expansion
`g_N(f)=g_N(0)+a_N f+O(f^2)`, with `a_N>0`, then near a loss of the
fold-adjoining stable component,

$$
f_c(N)\simeq -g_N(0)/a_N,\qquad
\Delta\epsilon(N)\ \propto\ [g_N(0)]^2
$$

on the side with `g_N(0)<0`. This last estimate is conditional: mode
switches can make the spectral abscissa nonsmooth, and no monotonicity in
the discrete chain length is assumed. A positive fold-limit growth
excludes a stable interval touching that fold by continuity; it does not
exclude detached stable intervals elsewhere on the branch.

Numerical size comparisons and their scope are recorded in
[the stability report](../docs/uniform_topological_segment_stability.md).
