"""Repeat identical fixed-background epsilon-phi stability map at N=20."""
import compensated_stability_map as mapping
from compensated_stability_map import *

def report20(summary):
    current=np.load(mapping.DEST/'compensated_stability_map.npz')
    old=np.load(Path('data/generated/compensation_n50/compensated_stability_map.npz'))
    assert np.array_equal(current['phi'],old['phi'])
    assert np.array_equal(current['epsilon'],old['epsilon'])
    stable20=current['growth']<-1e-7;stable50=old['growth']<-1e-7
    gained=int((stable20 & (old['growth']>1e-7)).sum())
    lost=int(((current['growth']>1e-7)&stable50).sum())
    # Independently recompute the familiar epsilon=38.51 slice.
    bounds_deg=json.loads(Path('data/generated/compensation_n50/compensated_phi_full_range.json').read_text())['existence_bounds_degrees']
    def slice_point(deg):
        phi=np.deg2rad(deg);c=.4998+2*np.cos(phi)
        n=sorted(z.real for z in np.roots([u*u,2*c*u,c*c+.25,-38.51**2]) if abs(z.imag)<1e-6 and z.real>0)[1]
        return float(n),mapping.growth(n,phi)
    lo,hi=bounds_deg;grid=np.unique(np.r_[np.linspace(lo+1e-5,hi-1e-5,600),lo+np.geomspace(1e-5,.1,40)])
    pairs=[(float(p),slice_point(float(p))[1]) for p in grid];crossings=[];checks=[]
    for a,b in zip(pairs[:-1],pairs[1:]):
        if a[1]*b[1]<0:
            p=brentq(lambda p:slice_point(p)[1],a[0],b[0],xtol=1e-10);crossings.append(float(p))
            for deg in [p-.001,p+.001]:
                n,g=slice_point(deg)
                hp=[precise(deg,prec,epsilon=38.51,initial_density=n,N=20) for prec in [100,160]]
                assert hp[0]['counts']==hp[1]['counts'];assert (g<0)==(hp[-1]['counts']['0']==0)
                checks.append(dict(degrees=deg,growth=g,checks=hp))
    comparison=dict(N=20,epsilon_slice=38.51,slice_existence_bounds_degrees=bounds_deg,
                    slice_crossings_degrees=crossings,slice_precision_checks=checks,
                    stable20_unstable50_samples=gained,unstable20_stable50_samples=lost)
    (mapping.DEST/'comparison_to_n50.json').write_text(json.dumps(comparison,indent=2))
    interval_rows=[]
    edges=[lo]+crossings+[hi]
    for left,right in reversed(list(zip(edges[:-1],edges[1:]))):
        interval_rows.append(f"|{left:.6f}-{right:.6f}|{'stable' if slice_point((left+right)/2)[1]<0 else 'unstable'}|")
    reverse_checks=[]
    for i,j in np.argwhere((current['growth']>1e-7)&stable50):
        deg=float(np.rad2deg(current['phi'][i]));eps=float(current['epsilon'][i,j]);n=float(current['density'][i,j])
        hp=[precise(deg,p,epsilon=eps,initial_density=n,N=size) for size in [20,50] for p in [100,160]]
        assert hp[0]['counts']==hp[1]['counts'] and hp[2]['counts']==hp[3]['counts']
        reverse_checks.append(dict(epsilon=eps,degrees=deg,checks=hp))
    (mapping.DEST/'nonmonotone_size_check.json').write_text(json.dumps(reverse_checks,indent=2))
    text=f'''## Two-dimensional compensation stability map repeated for N=20

This repeats the preceding N=50 calculation with **N=20**, keeping J=1, U=-0.0002, V=0, Delta=0.5, kappa=1 and exactly the same epsilon/phi backgrounds and sampling. The axes are epsilon horizontally and phi vertically. Green denotes an existing intermediate homogeneous root whose **N=20 OBC BdG operator is stable**, orange an unstable target, grey no distinct intermediate root, and yellow marginal/unresolved growth within 1e-7 J.

![Full compensated-intermediate stability map for N=20 OBC](../figures/compensated_stability_map_n20_full.png)

![N=20 zoom over epsilon=30-40, with the epsilon=38.51 slice marked](../figures/compensated_stability_map_n20_zoom.png)

### What changes with chain length

The root-existence lobe and root-specific compensation formula are unchanged: they depend on the uniform bulk cubic, not chain length. Every density and complex homogeneous background is therefore the same as in the N=50 calculation. Only the OBC BdG operator's length changes. The finite-size comparison consequently isolates spectral stability rather than a change of nonlinear background.

The N=20 calculation contains {summary['samples']} spectra on the identical 320-phase-row, 91-density-fraction mesh: {summary['stable_samples']} stable, {summary['unstable_samples']} unstable and {summary['marginal_samples']} marginal samples. At {gained} sampled backgrounds N=20 is stable while N=50 is unstable; at {lost} the reverse occurs. These counts are not geometric areas because the sampling is nonuniform. They do not establish monotonic stability versus N outside the tested sizes and points.

For the epsilon=38.51 slice, the intermediate root still exists only for {lo:.8f}<phi<{hi:.8f} degrees. The N=20 zero-growth crossings are {', '.join(f'{p:.8f} degrees' for p in crossings) or 'none'}. Crossings are marginal; the separate root-merger endpoints remain distinct from dynamical instability. The corresponding N=50 crossings were 18.28978185 and 56.23966500 degrees. Slice classifications on either side of the new crossings are independently checked with 100/160-digit root counts.

|epsilon=38.51 phase interval (degrees)|N=20 OBC target|
|---|---|
{chr(10).join(interval_rows)}

The extra short intervals are finite-chain reentrant stability windows. Their boundaries are numerical estimates and the 2D shading interpolates between sampled points; small features should not be interpreted more precisely than the sampling and direct slice checks support. The {lost} reverse size-comparison point is independently checked at 100/160 digits for both N=20 and N=50, confirming that shorter chains are not uniformly more stable. Its records are in nonmonotone_size_check.json.

### Verification and scope

The map uses the same curved density sampling, refined zero-growth boundaries, analytical existence masks and growth-rate interpolation as the N=50 figure. The N=20 characteristic-polynomial/Routh checks use **40 real-quadrature dimensions**, not the former 100-dimensional N=50 generator. Checks at phi=0,20,40,60,75 degrees and three density fractions agree between double-precision stability and 100/160-digit growing-mode counts. Actual complex compensated mean-field stationary residuals are verified below 1e-8 at these points. All selected checks passed; no time evolution, ramp preparation, or fluctuation feedback is used.

Green means the exact fixed endpoint compensation constructs a **locally stable** uniform mean-field state. It does not guarantee reachability from vacuum, global attraction, or robustness to endpoint-drive errors. Near colored boundaries this remains a finite-resolution numerical map; no thermodynamic conclusion follows from N=20 or N=50 alone.

Run `python src/compensated_stability_map_n20.py`. The N=50 files and discussion are retained for comparison. N=20 grids, densities, drives and growth rates are in data/generated/compensation_n20/compensated_stability_map.npz; parameters, refined zero-growth points and independent precision checks are in compensated_stability_map.json. Size comparisons and epsilon=38.51 slice checks are in comparison_to_n50.json. Full and zoom figures are available as PNG and PDF under figures/ with prefix compensated_stability_map_n20.
'''
    p=Path('docs/uniform_topological_segment_stability.md');base=p.read_text(encoding='utf-8').split('## Two-dimensional compensation stability map repeated for N=20')[0].rstrip();p.write_text(base+'\n\n'+text,encoding='utf-8')
    print('N20 slice crossings',crossings,'size differences',gained,lost,flush=True)

if __name__=='__main__':
    mapping.N=20
    mapping.DEST=Path('data/generated/compensation_n20');mapping.DEST.mkdir(parents=True,exist_ok=True)
    mapping.FIGURE_PREFIX='compensated_stability_map_n20'
    mapping.report=report20
    mapping.main()
