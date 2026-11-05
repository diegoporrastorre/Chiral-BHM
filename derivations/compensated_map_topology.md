# Zero-frequency topology of the compensated intermediate-state map

This extends the N=50, V=0, Delta=0.5, U=-0.0002, J=kappa=1
epsilon-phi map. It uses the same homogeneous stationary backgrounds and
fixed endpoint-drive compensation. No Gaussian feedback is included.

Use the SVD/extended chiral definition already documented in
`svd_topology_and_bulk_stability.md` (including its supplemental-source
provenance). At omega_ref=0, Q(k)=-H(k), the doubled Hermitian operator
is [[0,Q],[Q^dagger,0]], and its gap is the smallest singular value of Q.
Its chiral winding is

$$
\nu(0)=\frac{1}{2\pi i}\int_{-\pi}^{\pi} dk\,
\partial_k\log\det[-H(k)].
$$

The convention is forward hopping J exp(+i phi), Fourier exp(ikj), and
time evolution exp(-i omega t). A periodic bulk continuation uses the
same closing-bond phase, with total twist N phi. This is a bulk invariant
of the homogeneous coefficients; the stability plotted alongside it is
that of the finite compensated OBC chain.

## Exact result

Put a=Delta+U=0.4998, h=2 cos(phi), t=0.0004 n, gamma=0.5.
Then the effective detuning is d=a-2t and the pairing magnitude is g=t.

$$
\det[-H(k)]=(2\sin\phi\sin k-i\gamma)^2
-(d+2\cos\phi\cos k)^2+g^2.
$$

For 0<phi<pi/2, its imaginary part has a fixed negative sign for
0<k<pi and a positive sign for -pi<k<0. Define

$$
X_0=t^2-\gamma^2-(a+h-2t)^2,\qquad
X_\pi=t^2-\gamma^2-(a-h-2t)^2.
$$

Away from zero crossings, nu=(sign(X_pi)-sign(X_0))/2.
For the intermediate root, the cubic input-output function has negative
derivative with respect to n. Since X_0 is minus that derivative,
X_0>0 throughout the distinct intermediate branch.

Moreover, X_pi is strictly negative throughout the existence lobe.
To see this, set b=a-h. Then

$$
X_\pi=-3t^2+4bt-b^2-\gamma^2.
$$

If b<=0 this is negative for t>=0. If b>0, completing the square gives
X_pi=-3(t-2b/3)^2+b^2/3-gamma^2. Here 0<b<=a=0.4998,
so b^2/3-gamma^2<0. The existence lobe has phi<79.448881 degrees,
and hence h>0 as required.

Therefore **every distinct intermediate homogeneous state with phi>0
in this map has nu(0)=-1**, whether its finite OBC spectrum is stable or
unstable. The sign follows our spatial orientation; nonzero winding is
the topological conclusion.

At phi=0 the determinant is real and changes sign between k=0 and pi,
so it vanishes at an intermediate momentum: the point gap closes and
the invariant is undefined. At the root-merger curves X_0=0, the point
gap also closes. Neither boundary is included in the gapped conclusion.
Grey outside the existence lobe means no distinct intermediate root;
it says nothing about the topology of the other stationary branches.

## Numerical calculation and plot

`src/compensated_topology_map.py` evaluates the analytic invariant on the
original 320-by-91 density/phase mesh and independently sums determinant
phase increments on 4096 and 8192 momenta for every background. Both grids
include k=0 and k=-pi, resolving the real-axis crossings even close to
the root mergers. The analytic sign argument certifies the open point
gap; a finite momentum-grid minimum alone would not do so.

The original saved numerical grid was missing, so the script regenerates
the same stationary backgrounds and OBC spectra and refines their
zero-growth contours. The earlier report's selected high-precision
stability checks are not rerun or presented as new checks here.
Stability shading is finite-resolution interpolation, with marginal
growth tolerance 1e-7 J; no thermodynamic stability claim follows.

New PNG and vector SVG plots are saved under
`figures/compensated_stability_topology_{zoom,full}`. Computed arrays and
winding-check statistics are saved under
`data/generated/compensation_n50/compensated_topology_map.{npz,json}`.
The original figures are preserved.

Validation on 2026-11-05: both winding resolutions matched the exact
invariant at all 29,120 backgrounds (zero mismatches). Maximum phase
increments were 2.57918 and 2.13140 radians, respectively, below pi.
The regenerated stability counts reproduce the original report exactly:
12,875 stable and 16,245 unstable samples, all with nu=-1. These are
nonuniform sample counts, not parameter-space area fractions.
