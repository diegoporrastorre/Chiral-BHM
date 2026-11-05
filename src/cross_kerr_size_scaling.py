"""Size scaling of compensated intermediate branches at V=4U.

All rates are generator growth rates max Im(omega), in units J=1.
Decimal Routh counts resolve spectra when double precision is unreliable.
Generated JSON caches include the precision and every tested parameter.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import json
from pathlib import Path
from decimal import Decimal as D, localcontext
import numpy as np
from scipy.linalg import eigvals
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cross_kerr_stable_trivial_map import background
from cross_kerr_trivial_search import operator, symbol
from verify_cross_kerr_trivial import decimal_generator
from banded_stability_decimal import characteristic_banded, count_right
from mean_field import Chain

OUT = Path('data/generated/cross_kerr_trivial/size_scaling')
OUT.mkdir(parents=True, exist_ok=True)
U = -.0002
V = 4*U


def parameters(degrees, fraction):
    phi = np.deg2rad(degrees)
    n, epsilon = background(phi, fraction)
    return dict(ratio=4, V=V, delta=3., kappa=1., phi=float(phi),
                fraction=float(fraction), density=float(n), epsilon=float(epsilon))


def double_rate(N, degrees, fraction):
    r = parameters(degrees, fraction)
    return float(eigvals(operator(N, r['density'], r['phi'], 3, 1, V),
                         check_finite=False).imag.max())


def drive_width(degrees, fraction):
    """Exact cubic fold expansion, avoiding subtraction of nearby drives."""
    c = 3+U+2*np.cos(np.deg2rad(degrees))
    rad = np.sqrt(c*c-.75)
    numerator = -4*rad**3*fraction**2*(3-2*fraction)/(27*2*(U+V))
    return float(numerator/(parameters(degrees, 0)['epsilon']+
                            parameters(degrees, fraction)['epsilon']))


def decimal_check(N, degrees, fraction, precision=None, rate=False):
    """Count growing modes; optionally bracket the rightmost growth to 2e-6.

    Recomputes the background in Decimal; two precisions are checked at
    selected points separately. This is not an interval-arithmetic proof.
    """
    precision = precision or N+80
    key = f'N{N}_p{degrees:g}_f{fraction:.14g}_d{precision}_r{int(rate)}'
    path = OUT/(key+'.json')
    if path.exists():
        return json.loads(path.read_text())
    r = parameters(degrees, fraction)
    with localcontext() as ctx:
        ctx.prec = precision
        A, n, phi = decimal_generator(r, precision, N)
        coeff = characteristic_banded(A, 3)
        counts = {str(s): count_right(coeff, D(str(s)))
                  for s in [-1e-7, 0, 1e-7]}
        result = dict(N=N, degrees=degrees, fraction=fraction,
                      precision=precision, counts=counts)
        if rate:
            approx = double_rate(N, degrees, fraction)
            lo, hi = D(str(approx-1e-6)), D(str(approx+1e-6))
            if count_right(coeff, lo) == 0 or count_right(coeff, hi) != 0:
                lo, hi = D('-.6'), D('2')
                assert count_right(coeff, lo) > 0 and count_right(coeff, hi) == 0
                while hi-lo > D('0.000002'):
                    mid = (lo+hi)/2
                    if count_right(coeff, mid):
                        lo = mid
                    else:
                        hi = mid
            result.update(growth_lower=float(lo), growth_upper=float(hi),
                          double_growth=approx)
    path.write_text(json.dumps(result, indent=2))
    print('checked', key, result, flush=True)
    return result


def boundary(N, degrees):
    """First instability connected to the lower-density branch merger.

    The fold itself is not a distinct intermediate solution. Zero width
    denotes absence of the stable component adjoining that fold only.
    """
    fold = decimal_check(N, degrees, 0, rate=True)
    if fold['counts']['0']:
        return dict(N=N, degrees=degrees, fraction_lower=0., fraction_upper=0.,
                    drive_width_lower=0., drive_width_upper=0., fold=fold)
    # At the longest size continue the already verified boundary instead
    # of using an ill-conditioned dense spectrum as the initial guess.
    if N >= 180 and degrees == 130 and (OUT/'progress.json').exists():
        previous = [row['topological'] for row in json.loads((OUT/'progress.json').read_text())
                    if row['N'] < N]
        if previous:
            guess = previous[-1]['fraction_lower']
            lo, hi = max(0., guess-.005), min(.5, guess+.005)
            if (decimal_check(N, degrees, lo)['counts']['0'] == 0 and
                    decimal_check(N, degrees, hi)['counts']['0'] > 0):
                while hi-lo > 2e-5:
                    mid = (lo+hi)/2
                    if decimal_check(N, degrees, mid)['counts']['0']:
                        hi = mid
                    else:
                        lo = mid
                return dict(N=N, degrees=degrees, fraction_lower=lo, fraction_upper=hi,
                            drive_width_lower=drive_width(degrees, lo),
                            drive_width_upper=drive_width(degrees, hi), fold=fold)
    # Use dense eigenvalues to propose a narrow bracket, then independently
    # test its sides. At ill-conditioned sizes fall back to Routh bisection.
    coarse = brentq(lambda f: double_rate(N, degrees, f), 0, .5, xtol=1e-11)
    tolerance = 1e-8 if degrees == 35 else 1e-6
    lo, hi = max(0., coarse-tolerance), coarse+tolerance
    lc = decimal_check(N, degrees, lo)
    hc = decimal_check(N, degrees, hi)
    if lc['counts']['0'] or not hc['counts']['0']:
        lo, hi = 0., .5
        assert decimal_check(N, degrees, hi)['counts']['0'] > 0
        while hi-lo > 2*tolerance:
            mid = (lo+hi)/2
            if decimal_check(N, degrees, mid)['counts']['0']:
                hi = mid
            else:
                lo = mid
    return dict(N=N, degrees=degrees, fraction_lower=lo, fraction_upper=hi,
                drive_width_lower=drive_width(degrees, lo),
                drive_width_upper=drive_width(degrees, hi), fold=fold)


def validate_background(degrees, fraction):
    r = parameters(degrees, fraction)
    n, p, eps = r['density'], r['phi'], r['epsilon']
    alpha = -eps/(3+U+2*np.cos(p)+2*(U+V)*n-.5j)
    model = Chain(N=50, drive=eps, delta=3, kappa=1, U=U, V=V, phi=p)
    model.pump = model.pump.astype(complex)
    model.pump[0] += (np.exp(1j*p)+V*n)*alpha
    model.pump[-1] += (np.exp(-1j*p)+V*n)*alpha
    residual = float(abs(model.force(np.full(50, alpha))).max())
    actual = eigvals(model.bogoliubov(np.full(50, alpha))).imag.max()
    assert residual < 1e-10
    assert abs(actual-double_rate(50, degrees, fraction)) < 1e-8
    winding = []
    for nk in [8192, 16384]:
        k = np.linspace(-np.pi, np.pi, nk, endpoint=False)
        a, d, b = symbol(k, n, p, 3, 1, V)
        H = np.zeros((nk, 2, 2), complex)
        H[:, 0, 0], H[:, 1, 1] = a+d, a-d
        H[:, 0, 1], H[:, 1, 0] = b, -b
        det = np.linalg.det(-H)
        winding.append(dict(nk=nk, winding=float(np.angle(np.roll(det, -1)/det).sum()/(2*np.pi)),
                            gap=float(np.linalg.svd(H, compute_uv=False)[:, -1].min())))
    return dict(**r, degrees=degrees, residual=residual, winding=winding)


def main():
    sizes = [20, 35, 50, 65, 75, 80, 81, 82, 84, 90, 100, 128, 160, 200]
    rows = []
    for N in sizes:
        row = dict(N=N, trivial=boundary(N, 35), topological=boundary(N, 130),
                   fixed_trivial=decimal_check(N, 35, .001, rate=True),
                   fixed_topological=decimal_check(N, 130, .1, rate=True))
        rows.append(row)
        (OUT/'progress.json').write_text(json.dumps(rows, indent=2))
        print('completed size', N, flush=True)
    refined = []
    for N, degrees, fraction in [(50, 35, .001), (80, 35, 0), (81, 35, 0), (82, 35, 0),
                                  (84, 35, 0), (200, 35, .001),
                                  (50, 130, .1), (200, 130, .1)]:
        a = decimal_check(N, degrees, fraction, rate=True)
        b = decimal_check(N, degrees, fraction, precision=N+140, rate=True)
        assert a['counts'] == b['counts']
        assert max(a['growth_lower'], b['growth_lower']) <= min(a['growth_upper'], b['growth_upper'])
        refined.append(b)
    # Verify two separated interior backgrounds at the largest size:
    # a nonzero interval is supported independently of the boundary fit.
    for fraction in [.05, .15]:
        refined.append(decimal_check(200, 130, fraction, precision=340, rate=True))
    edge = rows[-1]['topological']
    left = decimal_check(200, 130, edge['fraction_lower'], precision=340)
    right = decimal_check(200, 130, edge['fraction_upper'], precision=340)
    assert left['counts']['0'] == 0 and right['counts']['0'] > 0
    refined.extend([left, right])
    # Sample the inferred connected intervals, including near their edges.
    # These numerical checks do not constitute an interval certificate.
    interior_checks = []
    for row in rows:
        if row['N'] not in [50, 80, 200]:
            continue
        for kind, degrees in [('trivial', 35), ('topological', 130)]:
            edge = row[kind]['fraction_lower']
            if not edge:
                continue
            for scale in [.1, .5, .9]:
                item = decimal_check(row['N'], degrees, edge*scale)
                assert item['counts']['0'] == 0
                interior_checks.append(item)
    result = dict(parameters=dict(J=1, U=U, V=V, delta=3, kappa=1,
                                   omega_ref=0, compensation='drive only'),
                  backgrounds=[validate_background(35, .001), validate_background(130, .1)],
                  rows=rows, refined=refined, interior_checks=interior_checks)
    (OUT/'summary.json').write_text(json.dumps(result, indent=2))
    plot(result)


def plot(result):
    rows = result['rows']; sizes = np.array([r['N'] for r in rows])
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.3))
    fig.subplots_adjust(left=.065, right=.985, bottom=.21, top=.80, wspace=.30)
    colors = dict(trivial='#a542a5', topological='#26866d')
    for kind, label in [('trivial', r'Trivial: $\phi=35^\circ$, $\nu(0)=0$'),
                        ('topological', r'Topological: $\phi=130^\circ$, $\nu(0)=-1$')]:
        color = colors[kind]
        lo = np.array([r[kind]['drive_width_lower'] for r in rows])
        hi = np.array([r[kind]['drive_width_upper'] for r in rows])
        width = (lo+hi)/2
        axes[0].plot(sizes, np.where(width>0, width, np.nan), 'o-', color=color, label=label)
        axes[0].fill_between(sizes, lo, hi, color=color, alpha=.25)
        for ax, source in [(axes[1], 'fixed_'+kind), (axes[2], kind)]:
            rates = [r[source] if ax is axes[1] else r[source]['fold'] for r in rows]
            gl = np.array([r['growth_lower'] for r in rates])
            gh = np.array([r['growth_upper'] for r in rates])
            ax.plot(sizes, (gl+gh)/2, 'o-', color=color)
            ax.fill_between(sizes, gl, gh, color=color, alpha=.25)
    axes[0].set(yscale='log', ylabel=r'Stable drive width $\Delta\epsilon/J$',
                title='Stable component adjoining the fold')
    axes[0].annotate('Trivial component absent\nat sampled '+r'$N\geq81$', xy=(81, 4e-7),
                     xytext=(105, 2e-5), arrowprops=dict(arrowstyle='->'), fontsize=9)
    axes[1].set(ylabel=r'Growth $\max\,\mathrm{Im}\,\omega/J$',
                title='Fixed distinct intermediate backgrounds')
    axes[1].text(.97, .06, r'Trivial $f=0.001$'+'\n'+r'Topological $f=0.1$',
                 transform=axes[1].transAxes, ha='right', fontsize=10)
    axes[2].set(ylabel=r'Fold-limit growth $\max\,\mathrm{Im}\,\omega/J$',
                title=r'Approaching the merger ($f\to0^+$)')
    for ax in axes:
        ax.set_xlabel('Chain length $N$'); ax.grid(alpha=.2)
    for ax in axes[1:]:
        ax.axhline(0, color='black', ls='--', lw=1)
    fig.suptitle(r'Compensated OBC size scaling: $J=\kappa=1$, $\Delta=3$, $U=-0.0002$, $V=4U$', fontsize=13)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(.5, .02),
               ncol=2, fontsize=10, frameon=False)
    for ext in ['png', 'svg']:
        fig.savefig(f'figures/cross_kerr_size_scaling.{ext}', dpi=200)


if __name__ == '__main__':
    main()
