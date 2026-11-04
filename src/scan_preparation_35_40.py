"""Complete N=50 preparation scan and replace superseded drive discussions."""
from scan_preparation_37_40 import *
from scipy.optimize import brentq

FILE35=DEST/'intermediate_preparation_35_40.json'

def save(r):
    import time
    tmp=FILE35.with_suffix('.writing');tmp.write_text(json.dumps(r))
    for attempt in range(25):
        try:
            tmp.replace(FILE35);return
        except PermissionError:
            if attempt==24:raise
            time.sleep(.2)

def main():
    result=json.loads(FILE35.read_text()) if FILE35.exists() else json.loads(FILE.read_text())
    pending=FILE35.with_suffix('.writing')
    if pending.exists():
        recovery=json.loads(pending.read_text())
        if len(recovery['records'])>len(result['records']):result=recovery
    grid=np.round(np.linspace(35,40,101),6).tolist()
    def add(values):
        done={r['epsilon'] for r in result['records']}
        for eps in sorted(set(values)-done):
            r=calculate(float(eps));result['records'].append(r);save(result);print(eps,outcome(r),flush=True)
    add(grid)
    coarse=sorted([r for r in result['records'] if r['epsilon'] in grid],key=lambda r:r['epsilon'])
    refine=[]
    for a,b in zip(coarse[:-1],coarse[1:]):
        if outcome(a)!=outcome(b):refine.extend(np.round(np.arange(a['epsilon']+.01,b['epsilon'],.01),6).tolist())
    add(refine)
    rr=sorted(result['records'],key=lambda r:r['epsilon']);selected={35.,38.51}
    for i in [0,1]:
        for a,b in zip(rr[:-1],rr[1:]):
            if a.get('three_roots') and b.get('three_roots') and outcome(a)[i]!=outcome(b)[i]:selected.update([a['epsilon'],b['epsilon']])
    done={c['epsilon'] for c in result['checks']}
    for eps in sorted(selected-done):
        r=next(r for r in rr if r['epsilon']==eps);a0=complex(*r['alpha0'])
        a,base=trajectory(Chain(N=50,drive=eps),np.zeros(50,complex),2400,2e-12,2e-13)
        check=dict(epsilon=eps,baseline_end_force=base['end_force'],runs=[])
        for T in [1000,10000]:check['runs'].append(compact(ramp(a.copy(),a0,T,2e-12,2e-13,epsilon=eps)))
        check['precision_checks']=[precise_count(50,.4998-.0008*r['density_target'],.0004*r['density_target'],p) for p in [100,160]]
        result['checks'].append(check);save(result);print('tight',eps,[z['hold']['target_reached'] for z in check['runs']],flush=True)
    # Exact turning points of the uniform bulk cubic (not a dynamical boundary).
    c=1.4998;u=-.0004
    turns=sorted([(float(np.sqrt(n*((c+u*n)**2+.25))),float(n)) for n in np.roots([3*u*u,4*c*u,c*c+.25])])
    result['bulk_turning_points']=[dict(epsilon=e,density=n) for e,n in turns]
    # Document intermediate uniform OBC spectra below the requested dynamics window.
    below=[]
    for eps in [30.17,31.,32.,33.,34.]:
        n=sorted(np.roots([u*u,2*c*u,c*c+.25,-eps*eps]).real)[1];a0=-eps/(c+u*n-.5j)
        z=diagnose(compensated(1,a0,eps),np.full(50,a0));below.append(dict(epsilon=eps,density_target=float(n),growth=z['growth'],winding=winding(1,np.pi/3,.4998-.0008*n,.0004*n,.5)['winding'],residual=z['residual']))
    result['below_window_target_only']=below;save(result);draw35(result);report35(result)

def draw35(result):
    rr=sorted([r for r in result['records'] if r.get('three_roots')],key=lambda r:r['epsilon'])
    fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for i,T in enumerate([1000,10000]):
        ax[0,0].semilogy([r['epsilon'] for r in rr],[r['runs'][i]['hold']['target_amplitude_error'] for r in rr],'.-',label=f'T={T}')
        ax[0,1].plot([r['epsilon'] for r in rr],[np.mean(r['runs'][i]['hold']['density']) for r in rr],'.-',label=f'T={T}')
    ax[0,0].axhline(1e-3,color='k',ls='--');ax[0,0].legend();ax[0,0].set(xlabel='epsilon / J',ylabel='max complex amplitude error')
    ax[0,1].plot([r['epsilon'] for r in rr],[r['density_target'] for r in rr],'k--',label='uniform intermediate target');ax[0,1].legend();ax[0,1].set(xlabel='epsilon / J',ylabel='final mean density')
    ax[1,0].plot([r['epsilon'] for r in rr],[r['target']['growth'] for r in rr],'.-');ax[1,0].axhline(0,color='k',ls='--');ax[1,0].set(xlabel='epsilon / J',ylabel='uniform target OBC growth / J')
    for eps in [35,35.5,36,36.5,37,38.51,39.85]:
        r=next(r for r in rr if r['epsilon']==eps);ax[1,1].plot(range(1,51),r['runs'][1]['hold']['density'],label=f'{eps:g}')
    ax[1,1].legend(title='epsilon / J',fontsize=8);ax[1,1].set(xlabel='site j',ylabel='final density, T=10000')
    fig.savefig('figures/intermediate_preparation_35_40.png',dpi=180)
    fig,ax=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    dense=np.linspace(30.16063463199232+.00001,39.870659871545925-.00001,800)
    densities=[sorted(np.roots([.0004**2,2*1.4998*(-.0004),1.4998**2+.25,-e*e]).real) for e in dense]
    for i,label in enumerate(['low','intermediate','high']):ax[0].plot(dense,np.array(densities)[:,i],label=label)
    ax[0].axvspan(35,40,color='.8',alpha=.3,label='requested scan');ax[0].set(xlabel='epsilon / J',ylabel='uniform bulk density');ax[0].legend()
    samples=result['below_window_target_only'];ax[1].plot([z['epsilon'] for z in samples],[z['growth'] for z in samples],'o-',label='target-only checks below 35');ax[1].plot([z['epsilon'] for z in rr],[z['target']['growth'] for z in rr],label='targets in preparation scan');ax[1].axhline(0,color='k',ls='--');ax[1].set(xlabel='epsilon / J',ylabel='uniform OBC N=50 growth / J');ax[1].legend(fontsize=8)
    fig.savefig('figures/intermediate_bulk_existence.png',dpi=180)

def report35(result):
    rr=sorted([r for r in result['records'] if r.get('three_roots')],key=lambda r:r['epsilon'])
    groups=[]
    for i,T in enumerate([1000,10000]):
        for lo,hi,success,count in ranges(rr,i):groups.append(f"|{T}|{lo:.2f}-{hi:.2f}|{count}|{'yes' if success else 'NO'}|")
    checks='; '.join(f"{c['epsilon']:.2f}: "+'/'.join('yes' if z['hold']['target_reached'] else 'NO' for z in c['runs']) for c in sorted(result['checks'],key=lambda c:c['epsilon']))
    details=[f"|{r['epsilon']:.2f}|{r['density_target']:.3f}|{r['target']['growth']:.6f}|{r['local_winding']['winding']}|"+'|'.join('yes' if z['hold']['target_reached'] else 'NO' for z in r['runs'])+'|' for r in rr]
    below=[f"|{r['epsilon']:.2f}|{r['density_target']:.3f}|{r['growth']:.6f}|{r['winding']}|" for r in result['below_window_target_only']]
    missing=sorted(r['epsilon'] for r in result['records'] if not r.get('three_roots'))
    lo,hi=result['bulk_turning_points']
    text=f'''## Complete drive scan: epsilon=35 to 40, N=50

This section replaces the earlier narrow and 37-40 drive-scan discussions. The three separate questions are **existence of the compensated uniform intermediate solution**, its **linear OBC stability and topology**, and **reachability with the specified vacuum-relaxation/ramp protocol**. A failed preparation does not imply nonexistence.

### Existence: the lower limit is below the requested interval

For the paper's U n^2 convention, the uniform bulk equation is

`epsilon^2 = n [(c+u n)^2+gamma^2]`,

where c=Delta+U+2J cos(phi)=1.4998, u=2U=-0.0004, gamma=kappa/2=0.5. Its turning points satisfy

`3u^2 n^2+4cu n+c^2+gamma^2=0`.

The resulting distinct three-root interval is

**{lo['epsilon']:.8f} < epsilon < {hi['epsilon']:.8f}.**

At the lower endpoint the intermediate and high roots merge at density {lo['density']:.6f}; at the upper endpoint the intermediate and low roots merge at density {hi['density']:.6f}. Exactly at a turning point two roots coincide, so there are not three distinct branches. Below the lower endpoint there is no separate intermediate uniform bulk root; above the upper endpoint there is no separate intermediate root either.

Root-specific endpoint compensation makes each such bulk root an exact uniform OBC stationary solution at any N. Therefore the **lower existence limit of the compensated uniform intermediate branch is epsilon={lo['epsilon']:.8f}, not a preparation threshold somewhere between 35 and 40**. This says nothing about disappearance of other inhomogeneous OBC stationary states.

![Uniform bulk roots, analytical turning points, and intermediate-target OBC stability checks](../figures/intermediate_bulk_existence.png)

To verify that stable topological targets can also exist below 35, we calculate their N=50 OBC spectra and winding at selected drives below the dynamical-scan window:

|epsilon / J|intermediate density|OBC growth / J|zero-frequency winding|
|---:|---:|---:|---:|
{chr(10).join(below)}

These are exact stationary-target and double-precision stability checks only; **no vacuum/ramp preparation was tested at these below-35 drives**. They do not prove stability at every drive down to the turning point, or at the turning point itself.

### Complete preparation protocol and sampling

For N=50, Delta=0.5, J=1, U=-0.0002, V=0, kappa=1, phi=pi/3, scan 35<=epsilon<=40 every 0.05, retain the earlier fine central points, and refine coarse brackets where either ramp outcome or root existence changes to 0.01 spacing. There are {len(result['records'])} sampled drives, including {len(rr)} with three roots. No intermediate target exists at the sampled drives {', '.join(f'{e:.2f}' for e in missing)}.

For each eligible drive independently: start from vacuum, evolve the uncompensated chain for 1200/J, then add that drive's intermediate-root compensation with s(t)=3(t/T)^2-2(t/T)^3, T=1000/J or 10000/J. Hold the final compensated drives for 1200/J. The endpoint additive pumps are epsilon+s J exp(+i phi) alpha0 and epsilon+s J exp(-i phi) alpha0; bulk pumps stay epsilon. Targets include both their density and complex phase. No root solve is inserted into the dynamical preparation.

Success requires a stationary raw endpoint (force <1e-3, distance from residual-accepted polished root <1e-3) and maximum complex amplitude distance from the uniform target <1e-3, including edge sites. Root acceptance is max|f|<1e-8, root tolerance 1e-10. DOP853 uses rtol=2e-10, atol=2e-11. The stability convention is exp(-i omega t); growth is max Im(omega), with classification tolerance 1e-7 J.

### Preparation outcomes

These spans summarize consecutive sampled drives with the same outcome; they are not proofs of continuum intervals or absence of finer nonmonotone windows.

|ramp time T (1/J)|sampled drive span|number of samples|uniform intermediate target reached?|
|---:|---|---:|---|
{chr(10).join(groups)}

![Complete preparation scan over epsilon=35-40: errors, means, target stability, and final spatial profiles](../figures/intermediate_preparation_35_40.png)

All eligible target growth rates in the requested window are negative: largest {max(r['target']['growth'] for r in rr):.6f} J. Target zero-frequency winding labels are {sorted(set(r['local_winding']['winding'] for r in rr),key=str)}. Thus sampled preparation failures in this window are not losses of local target stability. All uncompensated vacuum preparations pass the stationary check: {all(r['baseline']['reached_stationary'] for r in rr)}; their maximum raw endpoint force is {max(r['baseline']['end_force'] for r in rr):.3g}. For failed preparations, consult the recorded final profile, final-state growth, and tail variation to distinguish competing stationary attractors from unresolved dynamics.

### Independent checks and limitations

Repeat the complete protocol at epsilon=35, the existing reference drives, and both sides of each sampled preparation-outcome change with initial relaxation 2400/J, rtol=2e-12, atol=2e-13. Ramp and hold durations stay fixed. Tight-repeat outcomes (T=1000/T=10000): {checks}. These repeats and the corresponding 100/160-digit target Routh counts are recorded with the data. Double-precision spectra and precision convergence are not interval certificates. The method does not establish global uniqueness, all basins, or an asymptotic adiabatic limit.

Both ramp times succeed at every sampled drive from 35.00 through 39.61. Both fail at the sampled drives from 39.62 through 39.87 while the uniform target remains stable, reproducing the upper preparation limitation. Since both protocols still succeed at epsilon=35, **no lower preparation failure boundary is found in the requested interval**. The lower existence boundary epsilon=30.16063463 is analytical and independent of this preparation result; preparation below 35 remains untested.

### All eligible sampled drives

|epsilon / J|intermediate density|target OBC growth / J|zero-frequency winding|T=1000 target reached?|T=10000 target reached?|
|---:|---:|---:|---:|---|---|
{chr(10).join(details)}

Reproduce with `python src/scan_preparation_35_40.py`. Generated records are data/generated/compensation_n50/intermediate_preparation_35_40.json. Identical previously calculated drives are reused; every newly calculated drive starts independently from vacuum. Records include target phase and density, baseline convergence, sampled ramp trajectories, final complex amplitudes, residuals, force, tail variation, stability, winding, and independent repeats. The data preserve superseded scans for provenance; their separate drive-scan narratives have been replaced here. No covariance equations or fluctuation feedback are used.
'''
    p=Path('docs/uniform_topological_segment_stability.md');base=p.read_text(encoding='utf-8')
    markers=['## Drive dependence of intermediate-state preparation at N=50','## Extended preparation scan: epsilon=37 to 40, N=50','## Complete drive scan: epsilon=35 to 40, N=50']
    cuts=[base.index(marker) for marker in markers if marker in base];base=base[:min(cuts)].rstrip() if cuts else base.rstrip()
    p.write_text(base+'\n\n'+text,encoding='utf-8')

if __name__=='__main__':main()
