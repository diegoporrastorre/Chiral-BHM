"""Extend independently prepared intermediate-phase protocols to epsilon=37..40."""
from scan_intermediate_preparation import *

FILE=DEST/'intermediate_preparation_37_40.json'

def calculate(eps):
    c=1.4998;u=-.0004
    densities=sorted(x.real for x in np.roots([u*u,2*c*u,c*c+.25,-eps**2]) if abs(x.imag)<1e-7 and x.real>=0)
    if len(densities)!=3:return dict(epsilon=eps,three_roots=False)
    n=densities[1];a0=-eps/(c+u*n-.5j);m=compensated(1,a0,eps)
    a,base=trajectory(Chain(N=50,drive=eps),np.zeros(50,complex))
    base.pop('times');base.pop('density_trajectory')
    rec=dict(epsilon=eps,three_roots=True,density_target=float(n),alpha0=[float(a0.real),float(a0.imag)],
             target=diagnose(m,np.full(50,a0)),local_winding=winding(1,np.pi/3,.4998-.0008*n,.0004*n,.5),baseline=base,runs=[])
    for T in [1000,10000]:rec['runs'].append(compact(ramp(a.copy(),a0,T,epsilon=eps)))
    return rec

def outcome(r):
    return tuple(z['hold']['target_reached'] for z in r['runs']) if r.get('three_roots') else None

def main():
    result=json.loads(FILE.read_text()) if FILE.exists() else dict(N=50,records=json.loads((DEST/'intermediate_preparation_scan.json').read_text())['records'],checks=[])
    grid=np.round(np.linspace(37,40,61),6).tolist()
    def add(grid):
        done={r['epsilon'] for r in result['records']}
        for eps in grid:
            if eps in done:continue
            r=calculate(float(eps));result['records'].append(r);FILE.write_text(json.dumps(result));print(eps,outcome(r),flush=True)
    add(grid)
    # Resolve changes of either ramp outcome on a 0.01 grid within coarse brackets.
    coarse=sorted([r for r in result['records'] if r['epsilon'] in grid],key=lambda r:r['epsilon'])
    refine=[]
    for a,b in zip(coarse[:-1],coarse[1:]):
        if outcome(a)!=outcome(b):refine.extend(np.round(np.arange(a['epsilon']+.01,b['epsilon'],.01),6).tolist())
    add(refine)
    rr=sorted(result['records'],key=lambda r:r['epsilon'])
    selected={37.,40.,38.51}
    for i in [0,1]:
        for a,b in zip(rr[:-1],rr[1:]):
            if a.get('three_roots') and b.get('three_roots') and outcome(a)[i]!=outcome(b)[i]:selected.update([a['epsilon'],b['epsilon']])
    done={c['epsilon'] for c in result['checks']}
    for eps in sorted(selected-done):
        r=next(r for r in rr if r['epsilon']==eps)
        if not r.get('three_roots'):continue
        a0=complex(*r['alpha0']);a,base=trajectory(Chain(N=50,drive=eps),np.zeros(50,complex),2400,2e-12,2e-13)
        check=dict(epsilon=eps,baseline_end_force=base['end_force'],runs=[])
        for T in [1000,10000]:check['runs'].append(compact(ramp(a.copy(),a0,T,2e-12,2e-13,epsilon=eps)))
        check['precision_checks']=[precise_count(50,.4998-.0008*r['density_target'],.0004*r['density_target'],p) for p in [100,160]]
        result['checks'].append(check);FILE.write_text(json.dumps(result));print('tight',eps,[z['hold']['target_reached'] for z in check['runs']],flush=True)
    draw(result);write_report(result)

def ranges(rr,i):
    # List sampled changes, not inferred continuum intervals.
    runs=[]
    for r in rr:
        value=outcome(r)[i]
        if not runs or runs[-1][2]!=value:runs.append([r['epsilon'],r['epsilon'],value,1])
        else:runs[-1][1]=r['epsilon'];runs[-1][3]+=1
    return runs

def draw(result):
    rr=sorted([r for r in result['records'] if r.get('three_roots')],key=lambda r:r['epsilon'])
    fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for i,T in enumerate([1000,10000]):
        ax[0,0].semilogy([r['epsilon'] for r in rr],[z['hold']['target_amplitude_error'] for r in rr for z in [r['runs'][i]]],'.-',label=f'T={T}')
        ax[0,1].plot([r['epsilon'] for r in rr],[np.mean(r['runs'][i]['hold']['density']) for r in rr],'.-',label=f'T={T}')
    ax[0,0].axhline(1e-3,color='k',ls='--');ax[0,0].set(xlabel='epsilon / J',ylabel='maximum complex amplitude error');ax[0,0].legend()
    ax[0,1].plot([r['epsilon'] for r in rr],[r['density_target'] for r in rr],'k--',label='intermediate target');ax[0,1].set(xlabel='epsilon / J',ylabel='final mean density');ax[0,1].legend()
    ax[1,0].plot([r['epsilon'] for r in rr],[r['target']['growth'] for r in rr],'.-');ax[1,0].axhline(0,color='k',ls='--');ax[1,0].set(xlabel='epsilon / J',ylabel='uniform target OBC growth / J')
    for eps in [37.,37.5,38.,38.51,39.,39.5,39.85]:
        r=next(r for r in rr if r['epsilon']==eps);ax[1,1].plot(range(1,51),r['runs'][1]['hold']['density'],label=f'{eps:g}')
    ax[1,1].set(xlabel='site j',ylabel='final density, T=10000');ax[1,1].legend(title='epsilon / J',fontsize=8)
    fig.savefig('figures/intermediate_preparation_37_40.png',dpi=180)

def write_report(result):
    rr=sorted([r for r in result['records'] if r.get('three_roots')],key=lambda r:r['epsilon'])
    lines=[]
    for i,T in enumerate([1000,10000]):
        for lo,hi,success,count in ranges(rr,i):lines.append(f"|{T}|{lo:.2f}-{hi:.2f}|{count}|{'yes' if success else 'NO'}|")
    checks='; '.join(f"{c['epsilon']:.2f}: "+'/'.join('yes' if z['hold']['target_reached'] else 'NO' for z in c['runs']) for c in sorted(result['checks'],key=lambda c:c['epsilon']))
    detail=[]
    for r in rr:
        detail.append(f"|{r['epsilon']:.2f}|{r['density_target']:.2f}|{r['target']['growth']:.5f}|{r['local_winding']['winding']}|"+'|'.join('yes' if z['hold']['target_reached'] else 'NO' for z in r['runs'])+'|')
    missing=sorted(r['epsilon'] for r in result['records'] if not r.get('three_roots'))
    text=f'''## Extended preparation scan: epsilon=37 to 40, N=50

We extend the same independently vacuum-initialized protocol over 37<=epsilon<=40 on a coarse 0.05 grid, retaining the earlier fine central samples. Whenever either ramp outcome changes between adjacent coarse points, that bracket is sampled every 0.01. There are {len(result['records'])} sampled drives in total, of which {len(rr)} have three distinct positive bulk roots and therefore an intermediate target. Target density, complex phase, endpoint compensation, uniform-target BdG stability, and zero-frequency winding are recomputed independently at each eligible drive.

The intermediate root does not exist throughout the requested interval. Writing c=Delta+U+2J cos(phi)=1.4998 and u=2U=-0.0004, the turning points obey 3u^2 n^2+4cu n+c^2+(kappa/2)^2=0. The low and intermediate roots merge at n=1479.250071 and epsilon=39.87065987. Above this upper saddle-node there is only one positive homogeneous bulk root, so the intermediate-root compensation is undefined. Sampled drives with no three-root target: {', '.join(f'{e:.2f}' for e in missing)}. This is an analytical bulk-root existence boundary, separate from preparation success and from OBC BdG stability of the compensated target below it.

Each run begins from vacuum with **no compensation**, relaxes for 1200/J, applies intermediate-root endpoint compensation with s=3(t/T)^2-2(t/T)^3, then holds the final drives fixed for 1200/J. Both T=1000/J and T=10000/J are tested. Other parameters remain N=50, Delta=0.5, J=1, U=-0.0002, V=0, kappa=1, phi=pi/3. Success requires stationarity and maximum complex amplitude error below 1e-3, including endpoints. These are finite-time protocol outcomes, not a proof of an adiabatic limit.

### Sampled outcome groups

The endpoints below describe consecutive **sampled** drives with the same outcome. Gaps between samples have not been exhaustively tested, and nonmonotone success between samples is not excluded.

|ramp time T (1/J)|sampled drive span|number of samples|uniform intermediate target reached?|
|---:|---|---:|---|
{chr(10).join(lines)}

![Preparation outcomes and final profiles over epsilon=37-40](../figures/intermediate_preparation_37_40.png)

The largest target growth is {max(r['target']['growth'] for r in rr):.6f} J; zero-frequency winding labels sampled are {sorted(set(r['local_winding']['winding'] for r in rr),key=str)}. This distinguishes a target instability from a finite-ramp failure to enter its basin. A failed preparation does not establish that the uniform target is absent or globally inaccessible. The largest uncompensated preparation endpoint force is {max(r['baseline']['end_force'] for r in rr):.3g}; all preparation endpoints pass the stationary test: {all(r['baseline']['reached_stationary'] for r in rr)}.

### Tolerance and stability checks

At interval endpoints, epsilon=38.51, and both sides of every sampled change of outcome, repeat the full preparation with initial relaxation 2400/J, rtol=2e-12 and atol=2e-13. Ramp times and final hold remain unchanged. Tight-repeat results (T=1000/T=10000): {checks}. Any mismatch with the base outcomes means that threshold should be regarded as numerically sensitive rather than established.

For the same selected target operators, 100- and 160-digit characteristic-polynomial/Routh counts are stored to check the number of growing eigenmodes. Double-precision spectra elsewhere are not interval certificates. In this scan, stable/unstable classification uses tolerance 1e-7 J; near-marginal targets require additional verification before drawing conclusions.

### All sampled drives

|epsilon / J|intermediate target density|target OBC growth / J|zero-frequency winding|T=1000 target reached?|T=10000 target reached?|
|---:|---:|---:|---:|---|---|
{chr(10).join(detail)}

Driver: `python src/scan_preparation_37_40.py`. Generated records are data/generated/compensation_n50/intermediate_preparation_37_40.json. Previously computed central drives are reused because their parameters, preparation, and tolerances are identical. New drives are calculated independently from vacuum; no drive continuation is used for preparation. Records include sampled ramp profiles, final complex amplitudes, target and approached-state stability, stationary residuals, raw endpoint forces, final temporal variation, root distances, and precision repeats. The final finite-ramp outcome is checked after integration; root polishing never replaces an initial state or steers the dynamics. The earlier narrower scan remains above for comparison.
'''
    p=Path('docs/uniform_topological_segment_stability.md');base=p.read_text(encoding='utf-8').split('## Extended preparation scan: epsilon=37 to 40, N=50')[0].rstrip();p.write_text(base+'\n\n'+text,encoding='utf-8')

if __name__=='__main__':main()
