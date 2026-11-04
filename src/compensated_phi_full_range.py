"""Full 0..pi intermediate-root existence and compensated OBC BdG stability."""
from compensated_phi_stability import point,precise,lower_turn
from compensation_n50_vacuum import *
from scipy.optimize import brentq

TARGET=DEST/'compensated_phi_full_range.json'

def upper_turn(degrees):
    c=.4998+2*np.cos(np.deg2rad(degrees));u=-.0004
    n=min(np.roots([3*u*u,4*c*u,c*c+.25]))
    return float(np.sqrt(n*((c+u*n)**2+.25)))

def main():
    lo=brentq(lambda p:lower_turn(p)-38.51,0,60,xtol=1e-12)
    hi=brentq(lambda p:upper_turn(p)-38.51,60,70,xtol=1e-12)
    grid=sorted(set(np.r_[np.linspace(0,180,1801),lo+np.geomspace(1e-5,.1,60),hi-np.geomspace(1e-5,.1,60)]),reverse=True)
    rows=[point(float(p)) for p in grid]
    valid=sorted([r for r in rows if r['intermediate_exists']],key=lambda r:r['degrees'])
    crossings=[]
    for a,b in zip(valid[:-1],valid[1:]):
        if a['growth']*b['growth']<0:
            root=brentq(lambda p:point(p)['growth'],a['degrees'],b['degrees'],xtol=1e-10)
            crossings.append(dict(degrees=float(root),phi=float(np.deg2rad(root)),scan_bracket=[a['degrees'],b['degrees']]))
    result=dict(N=50,epsilon=38.51,parameters=dict(J=1,U=-.0002,V=0,delta=.5,kappa=1),
                existence_bounds_degrees=[float(lo),float(hi)],existence_bounds_phi=[float(np.deg2rad(lo)),float(np.deg2rad(hi))],
                rows=rows,transitions=crossings,precision_checks=[])
    TARGET.write_text(json.dumps(result))
    selected={60,61,hi-.001,lo+.01,20,40,55}
    for c in crossings:selected.update([c['degrees']-.001,c['degrees']+.001])
    for degrees in sorted(selected,reverse=True):
        for precision in [100,160]:
            check=precise(degrees,precision);result['precision_checks'].append(check);TARGET.write_text(json.dumps(result));print(check,flush=True)
    draw(result);report(result)
    print('existence',lo,hi,'crossings',crossings,flush=True)

def draw(result):
    valid=sorted([r for r in result['rows'] if r['intermediate_exists']],key=lambda r:r['degrees'])
    lo,hi=result['existence_bounds_degrees'];deg=[r['degrees'] for r in valid]
    fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    ax[0,0].plot(deg,[r['density'] for r in valid]);ax[0,0].set(ylabel='intermediate density',title='Full requested phase interval')
    ax[0,1].plot(deg,[r['growth'] for r in valid]);ax[0,1].set(ylabel='OBC growth / J',title='Full requested phase interval')
    ax[1,0].plot(deg,[r['growth'] for r in valid]);ax[1,0].set(ylabel='OBC growth / J',title='Zoom on the existing intermediate branch')
    ax[1,1].semilogy(deg,[max(r['condition'],1) for r in valid]);ax[1,1].set(ylabel='dominant eigenvalue condition',title='Conditioning on the existing branch')
    for i,axes in enumerate(ax.flat):
        axes.axvspan(0,lo,color='.8',alpha=.5,label='no intermediate root');axes.axvspan(hi,180,color='.8',alpha=.5)
        axes.set(xlabel='phi (degrees)',xlim=(0,180) if i<2 else (17,63));axes.grid(alpha=.15)
        if i in [1,2]:
            axes.axhline(0,color='k',ls='--')
            for c in result['transitions']:axes.axvline(c['degrees'],ls=':',color='#b44c12')
    ax[0,0].legend(fontsize=8);fig.savefig('figures/compensated_phi_full_range.png',dpi=180)
    fig,axes=plt.subplots(1,3,figsize=(14,4),layout='constrained')
    for axes_one,limits in zip(axes,[(lo-.005,lo+.18),(55.8,56.7),(hi-.15,hi+.015)]):
        phases=np.linspace(max(limits[0],lo+1e-5),min(limits[1],hi-1e-5),140)
        data=[point(float(p)) for p in phases]
        axes_one.plot(phases,[r['growth'] for r in data]);axes_one.axhline(0,color='k',ls='--')
        axes_one.axvspan(0,lo,color='.8',alpha=.5);axes_one.axvspan(hi,180,color='.8',alpha=.5)
        for c in result['transitions']:axes_one.axvline(c['degrees'],ls=':',color='#b44c12')
        axes_one.set(xlim=limits,xlabel='phi (degrees)',ylabel='OBC growth / J');axes_one.grid(alpha=.15)
    fig.savefig('figures/compensated_phi_full_boundaries.png',dpi=180)

def report(result):
    lo,hi=result['existence_bounds_degrees'];lop,hip=result['existence_bounds_phi'];low_cross,high_cross=[r['degrees'] for r in result['transitions']]
    valid=[r for r in result['rows'] if r['intermediate_exists']]
    table=[]
    for degrees in [180,150,120,90,75,65,62,61.3,61,60,57,56,50,40,30,20,18.25,18,0]:
        r=point(degrees)
        table.append(f"|{degrees:g}|{r['density']:.3f}|{r['growth']:.6f}|{'stable' if r['growth']<-1e-7 else 'unstable' if r['growth']>1e-7 else 'marginal'}|" if r['intermediate_exists'] else f'|{degrees:g}|--|--|no intermediate root|')
    checks=[f"|{r['degrees']:.6f}|{r['counts']['0']}|" for r in result['precision_checks'] if r['precision']==160]
    text=f'''## Compensated intermediate-state stability over 0 <= phi <= pi

This complete phase description replaces the earlier pi/3-to-zero discussion. Fix **epsilon=38.51 and N=50**, J=1, U=-0.0002, V=0, Delta=0.5, kappa=1. Scan phi from pi (180 degrees) down to zero. This is a stationary mean-field and OBC Bogoliubov calculation only: **no adiabatic ramp or time evolution is simulated**.

### Existence comes before stability

At each phi recompute the positive real roots of

`epsilon^2 = n [(c(phi)+2U n)^2+(kappa/2)^2]`,

`c(phi)=Delta+U+2J cos(phi)=0.4998+2 cos(phi)`.

The intermediate root means the middle density among **three positive real roots**. It exists only in

**{lo:.8f} < phi < {hi:.8f} degrees**, or **{lop:.10f} < phi < {hip:.10f} radians**.

At the higher-phase endpoint {hi:.8f} degrees, the low- and intermediate-density roots merge. Thus, beginning at phi=pi and decreasing, there is initially only one positive bulk root, and a separate intermediate branch appears only below this endpoint. At the lower-phase endpoint {lo:.8f} degrees, the intermediate and high roots merge; below it only one positive root remains again. At each merger two roots coincide, so the endpoint is not part of the distinct three-root interval.

These existence bounds are solved from the two bulk cubic turning points, not inferred from a failed numerical steady-state solver. With u=2U, they satisfy 3u^2 n^2+4c(phi)u n+c(phi)^2+(kappa/2)^2=0 together with the stationary cubic. Positivity and reality are checked explicitly on the full phase grid. For phi where c<=sqrt(3) kappa/2, positive-density bistability is impossible; in particular no intermediate branch reappears at larger phi up to pi. The remaining single-root solution is not absent and is not assessed here as an intermediate phase.

### Compensation and Bogoliubov conventions

Where the intermediate root exists, set alpha0=-epsilon/[c(phi)+2Un-i kappa/2] and recompute the fixed endpoint additive drives:

`pump_1=epsilon+J exp(+i phi) alpha0`,

`pump_N=epsilon+J exp(-i phi) alpha0`.

Interior drive stays epsilon. This makes the target exactly uniform, including both edge sites; the maximum stationary-force residual on the full scan is {max(r['residual'] for r in valid):.3g}. The paper's U n^2 convention, including the +U shift, is retained. The density, complex phase, pairing, and endpoint compensation vary with phi; this is not a frozen-background quadratic sweep.

The fixed drives add no BdG derivatives. The OBC matrix is L=[[A,B],[-B*,-A*]], with A=T+diag(Delta+U+4Un-i kappa/2), B=2U alpha0^2 I. Directed nearest-neighbour hopping is T[j+1,j]=J exp(+i phi), sites 1..50, and no closing bond. Modes evolve as exp(-i omega t); the maximum Im(omega) is the growth rate. Negative is stable, positive unstable, and magnitudes below 1e-7 J are marginal/unresolved.

### Complete result, read from pi downward

|phase interval (degrees)|compensated intermediate state|
|---|---|
|{hi:.4f} < phi <= 180|no separate intermediate root|
|{high_cross:.4f} < phi < {hi:.4f}|stable|
|{low_cross:.4f} < phi < {high_cross:.4f}|unstable|
|{lo:.4f} < phi < {low_cross:.4f}|stable again, narrow interval|
|0 <= phi < {lo:.4f}|no separate intermediate root|

Zero-growth crossings occur at {high_cross:.8f} degrees and {low_cross:.8f} degrees and are marginal. The compensated intermediate state therefore does **not** exist near phi=pi or phi=0, and existence must not be confused with stability. Where it exists, its OBC BdG operator changes from stable to unstable and back to stable as phi decreases. These are finite-N=50 stability intervals, not thermodynamic conclusions or claims about the other roots.

![Full 0-to-pi scan: root existence, density, OBC stability and conditioning](../figures/compensated_phi_full_range.png)

![Resolved stability crossings and both branch-existence boundaries](../figures/compensated_phi_full_boundaries.png)

|phi (degrees)|intermediate density|OBC growth / J|classification|
|---:|---:|---:|---|
{chr(10).join(table)}

### Numerical resolution and independent precision checks

The full grid has 0.1-degree spacing over 0-180 degrees. Extra phases approach both existence endpoints from inside by logarithmically spaced distances 1e-5 to 0.1 degrees, resolving the narrow low-phase stable interval. Spectral-growth sign changes are refined with a scalar zero finder. Full complex OBC spectra, dominant left/right eigenvector condition numbers, residuals and sensitivity indicators are saved. Strong nonnormality means last digits of double-precision eigenvalues need care; stability-crossing digits are numerical estimates rather than certified bounds.

For representative phases, a point 0.001 degrees below the high-phase root merger, a point 0.01 degrees above the low-phase merger, and both sides of each stability crossing offset by 0.001 degrees, the background density and trigonometric coefficients are recomputed at high precision. The real-quadrature generator's banded characteristic polynomial and Routh recursion count growing roots at 100 and 160 digits, also checking shifts +/-1e-7. The 160-digit zero-shift counts are:

|phi (degrees)|growing generator modes|
|---:|---:|
{chr(10).join(checks)}

Counts at both precisions are saved. They independently check stable/unstable labels at selected phases, not every continuum point, and are precision-convergence evidence rather than interval-arithmetic certificates. No stationary root is extrapolated beyond a branch merger. No preparation or basin-of-attraction conclusion is inferred from this stability-only calculation.

### Reproduction and records

Run `python src/compensated_phi_full_range.py`. The complete dataset is data/generated/compensation_n50/compensated_phi_full_range.json, with all positive real bulk roots, existence flags across 0..pi, the intermediate complex amplitude and both complex endpoint pumps where defined, full OBC spectra, residuals, growth, conditioning, refined crossings, analytical merger locations and precision counts. The earlier partial-range dataset is retained for provenance; its discussion has been replaced by this complete description. No Gaussian covariance evolution or fluctuation feedback is used.
'''
    p=Path('docs/uniform_topological_segment_stability.md');base=p.read_text(encoding='utf-8')
    markers=['## Compensated intermediate-state stability versus phi at epsilon=38.51','## Compensated intermediate-state stability over 0 <= phi <= pi']
    cuts=[base.index(marker) for marker in markers if marker in base];base=base[:min(cuts)].rstrip() if cuts else base.rstrip()
    p.write_text(base+'\n\n'+text,encoding='utf-8')

if __name__=='__main__':main()
