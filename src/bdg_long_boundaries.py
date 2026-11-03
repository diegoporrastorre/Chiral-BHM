"""Long-chain scaling of the prescribed Rassaert-inspired BdG taper."""
from pathlib import Path
from decimal import Decimal, localcontext
import json, math, time
from banded_stability_decimal import threshold, envelope_generator, characteristic_banded, count_right

ROOT=Path(__file__).resolve().parents[1]
PATH=ROOT/'data/generated/bdg_long_boundaries.json'


def verified_threshold(n,width):
    precision=n+80
    start=time.perf_counter()
    # Locate a real leading-root candidate cheaply, then count ALL roots to
    # certify spectral ordering. Fall back to Routh bisection if necessary.
    with localcontext() as ctx:
        ctx.prec=precision
        coeff=characteristic_banded(envelope_generator(n,width,'paired'))
        def value(x):
            y=Decimal(0)
            for c in coeff[::-1]: y=y*x+c
            return y
        hi=Decimal(4); lo=hi
        step=Decimal('.01')
        while lo>step and value(lo)>=0:
            hi=lo; lo-=step
        if value(lo)<0:
            while hi-lo>Decimal('1e-5'):
                mid=(lo+hi)/2
                if value(mid)>0: hi=mid
                else: lo=mid
            counts=[count_right(coeff,lo),count_right(coeff,hi)]
        else: counts=[0,1]
        if counts[0]>0 and counts[1]==0:
            r=dict(N=n,width=width,kind='paired',precision=precision,
                   growth_lower=float(lo),growth_upper=float(hi),critical_kappa=float(lo+hi),
                   method='real-root candidate plus full Routh counts')
        else:
            r=threshold(n,width,'paired',precision,tolerance='1e-5')
            r['method']='full Routh bisection'
    with localcontext() as ctx:
        ctx.prec=precision+100
        coeff=characteristic_banded(envelope_generator(n,width,'paired'))
        lo=Decimal(str(r['growth_lower']))-Decimal('1e-10')
        hi=Decimal(str(r['growth_upper']))+Decimal('1e-10')
        counts=[count_right(coeff,lo),count_right(coeff,hi)]
        assert counts[0]>0 and counts[1]==0,(n,width,r,counts)
    r.update(check_precision=precision+100,check_counts=counts,
             kappa_lower=2*r['growth_lower'],kappa_upper=2*r['growth_upper'],
             center_gs=math.tanh(n/(2*width))**2,
             bulk_like_sites=sum(math.tanh(min(j,n+1-j)/width)**2>.99 for j in range(1,n+1)),
             elapsed_seconds=time.perf_counter()-start)
    return r


def run():
    data=json.loads(PATH.read_text()) if PATH.exists() else dict(
        parameters=dict(gs=1,gc=1,J=1,Delta=.1,phi='pi/2',P=0),
        profile='s=tanh(distance/N0), Delta=.1*s^2, gs=s^2, gc=s_j*s_(j+1); uniform kappa',rows=[])
    for n in [128,256,512]:
        for label,width in [('fixed8',8),('fraction1_8',n/8),('fraction1_5',n/5)]:
            if any(r['N']==n and r['series']==label for r in data['rows']): continue
            print('START',n,label,width,flush=True)
            r=verified_threshold(n,width); r['series']=label
            data['rows'].append(r)
            PATH.write_text(json.dumps(data,indent=2))
            print('DONE',n,label,r['critical_kappa'],r['check_counts'],r['elapsed_seconds'],flush=True)
    return data


def plot(data):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    labels={'fixed8':r'Fixed $N_0=8$','fraction1_8':r'$N_0=N/8$','fraction1_5':r'$N_0=N/5$'}
    fig,axs=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    for series,label in labels.items():
        rr=sorted([r for r in data['rows'] if r['series']==series],key=lambda r:r['N'])
        n=[r['N'] for r in rr]; kc=[r['critical_kappa'] for r in rr]
        axs[0].plot(n,kc,'o-',label=label)
        axs[1].plot(n,[(r['critical_kappa']-3.5)/2 for r in rr],'o-',label=label)
    axs[0].axhline(3.5,color='gray',ls=':',label=r'Example loss $\kappa=3.5J$')
    axs[0].axhline(2,color='k',ls='--',label='Original local bulk threshold')
    axs[0].axhline(2*np.sqrt(4.99),color='gray',ls='--',label='Sharp-OBC large-N threshold')
    axs[0].set(xlabel='Number of sites N',ylabel=r'Critical loss $\kappa_c/J$',title='Strict stability above each curve')
    axs[0].legend(fontsize=8)
    axs[1].axhline(0,color='k',lw=.8)
    axs[1].set(xlabel='Number of sites N',ylabel=r'Growth rate at $\kappa=3.5J$',title='Negative growth = stable')
    axs[1].legend(fontsize=8)
    fig.suptitle(r'Long chains: $\Delta=0.1J$, $g_s=g_c=J=1$, $\phi=\pi/2$; joint tanh taper')
    fig.savefig(ROOT/'figures/bdg_long_boundaries.png',dpi=180)
    plt.close(fig)


def write_report(data):
    lines=[r'''# Longer chains with Rassaert-inspired boundary profiles

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
''']
    lines.append('| N | Fixed N0=8 | N0=N/8 | N0=N/5 |\n|---:|---:|---:|---:|')
    for n in [128,256,512]:
        rr=[next(r for r in data['rows'] if r['N']==n and r['series']==s) for s in ['fixed8','fraction1_8','fraction1_5']]
        lines.append('| '+str(n)+' | '+' | '.join(f"{r['critical_kappa']:.4f}" for r in rr)+' |')
    lines.append(r'''
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
''')
    (ROOT/'derivations/bdg_long_boundary_scaling.md').write_text('\n'.join(lines))


if __name__=='__main__':
    result=run()
    plot(result)
    write_report(result)
