"""Complete low-drive ramp tests and consolidate epsilon=30..40 results."""
from scan_preparation_35_40 import *

LOW=DEST/'preparation_30_35'
LOW.mkdir(exist_ok=True)

def filename(eps):return LOW/f'epsilon_{eps:.6f}.json'

def store(path,value):
    import time
    temp=path.with_suffix('.writing');temp.write_text(json.dumps(value))
    for attempt in range(30):
        try:temp.replace(path);return
        except PermissionError:
            if attempt==29:raise
            time.sleep(.2)

def main():
    grid=np.round(np.linspace(30,35,101),6).tolist()
    # First broad pass, then fill the 0.05 grid. Every trajectory is independent.
    pilot=np.round(np.linspace(30,35,21),6).tolist()
    near=[30.16,30.1607,30.161,30.162,30.165,30.17,30.18,30.19]
    previous=json.loads(FILE35.read_text())
    for values in [pilot,grid,near]:
        for eps in values:
            p=filename(eps)
            if p.exists():continue
            known=next((r for r in previous['records'] if r['epsilon']==eps),None)
            r=known if known is not None else calculate(float(eps))
            store(p,r);print(eps,outcome(r),flush=True)
    coarse=[json.loads(filename(e).read_text()) for e in grid]
    extra=[]
    for a,b in zip(coarse[:-1],coarse[1:]):
        if outcome(a)!=outcome(b):extra.extend(np.round(np.arange(a['epsilon']+.01,b['epsilon'],.01),6).tolist())
    for eps in sorted(set(extra)):
        if filename(eps).exists():continue
        r=calculate(eps);store(filename(eps),r);print('refine',eps,outcome(r),flush=True)
    low=sorted([json.loads(p.read_text()) for p in LOW.glob('epsilon_*.json')],key=lambda r:r['epsilon'])
    selected={30.1607,30.17,31.,32.,33.,34.}
    for i in [0,1]:
        for a,b in zip(low[:-1],low[1:]):
            if a.get('three_roots') and b.get('three_roots') and outcome(a)[i]!=outcome(b)[i]:selected.update([a['epsilon'],b['epsilon']])
    checks=[]
    for eps in sorted(selected):
        p=LOW/f'check_{eps:.6f}.json'
        if p.exists():checks.append(json.loads(p.read_text()));continue
        r=next(r for r in low if r['epsilon']==eps);a0=complex(*r['alpha0'])
        a,base=trajectory(Chain(N=50,drive=eps),np.zeros(50,complex),2400,2e-12,2e-13)
        check=dict(epsilon=eps,baseline_end_force=base['end_force'],runs=[])
        for T in [1000,10000]:check['runs'].append(compact(ramp(a.copy(),a0,T,2e-12,2e-13,epsilon=eps)))
        check['precision_checks']=[precise_count(50,.4998-.0008*r['density_target'],.0004*r['density_target'],precision) for precision in [100,160]]
        checks.append(check);store(p,check);print('tight',eps,[z['hold']['target_reached'] for z in check['runs']],flush=True)
    mapping={r['epsilon']:r for r in previous['records']};mapping.update({r['epsilon']:r for r in low})
    result=dict(N=50,records=sorted(mapping.values(),key=lambda r:r['epsilon']),checks=previous['checks']+checks,
                bulk_turning_points=previous['bulk_turning_points'])
    # Keep large trajectory records in their source files; a compact index identifies them.
    index=dict(N=50,low_files=[str(p) for p in LOW.glob('epsilon_*.json')],upper_file=str(FILE35),
               low_check_files=[str(p) for p in LOW.glob('check_*.json')])
    store(DEST/'intermediate_preparation_30_40_index.json',index)
    draw30(result);report30(result)

def draw30(result):
    rr=[r for r in result['records'] if r.get('three_roots')]
    fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for i,T in enumerate([1000,10000]):
        ax[0,0].semilogy([r['epsilon'] for r in rr],[r['runs'][i]['hold']['target_amplitude_error'] for r in rr],'.-',label=f'T={T}')
        ax[0,1].plot([r['epsilon'] for r in rr],[np.mean(r['runs'][i]['hold']['density']) for r in rr],'.-',label=f'T={T}')
    ax[0,0].axhline(1e-3,color='k',ls='--');ax[0,0].legend();ax[0,0].set(xlabel='epsilon / J',ylabel='max complex amplitude error')
    ax[0,1].plot([r['epsilon'] for r in rr],[r['density_target'] for r in rr],'k--',label='intermediate target');ax[0,1].legend();ax[0,1].set(xlabel='epsilon / J',ylabel='final mean density')
    ax[1,0].plot([r['epsilon'] for r in rr],[r['target']['growth'] for r in rr],'.-');ax[1,0].axhline(0,color='k',ls='--');ax[1,0].set(xlabel='epsilon / J',ylabel='target OBC growth / J')
    for eps in [30.1607,30.17,30.5,31,32,33,34]:
        r=next(r for r in rr if r['epsilon']==eps);ax[1,1].plot(range(1,51),r['runs'][1]['hold']['density'],label=f'{eps:g}')
    ax[1,1].legend(title='epsilon / J',fontsize=8);ax[1,1].set(xlabel='site j',ylabel='final density, T=10000')
    fig.savefig('figures/intermediate_preparation_30_40.png',dpi=180)
    fig,ax=plt.subplots(1,2,figsize=(12,4),layout='constrained')
    low=[r for r in rr if r['epsilon']<=35]
    for i,T in enumerate([1000,10000]):
        ax[0].semilogy([r['epsilon'] for r in low],[r['runs'][i]['hold']['target_amplitude_error'] for r in low],'.-',label=f'T={T}')
        near=[r for r in low if r['epsilon']<=30.5]
        ax[1].semilogy([r['epsilon'] for r in near],[r['runs'][i]['hold']['target_amplitude_error'] for r in near],'.-',label=f'T={T}')
    for axes in ax:
        axes.axhline(1e-3,color='k',ls='--');axes.legend();axes.set(xlabel='epsilon / J',ylabel='max complex amplitude error')
    ax[0].set_title('Low-drive preparation');ax[1].set_title('Near lower root-existence boundary')
    fig.savefig('figures/intermediate_preparation_30_35_detail.png',dpi=180)

def report30(result):
    rr=[r for r in result['records'] if r.get('three_roots')];low=[r for r in rr if r['epsilon']<35]
    group=[]
    for i,T in enumerate([1000,10000]):
        for lo,hi,success,count in ranges(rr,i):group.append(f"|{T}|{lo:.6f}-{hi:.6f}|{count}|{'yes' if success else 'NO'}|")
    detail=[f"|{r['epsilon']:.6f}|{r['density_target']:.3f}|{r['target']['growth']:.6f}|{r['local_winding']['winding']}|"+'|'.join('yes' if z['hold']['target_reached'] else 'NO' for z in r['runs'])+'|' for r in rr]
    checks='; '.join(f"{c['epsilon']:.6f}: "+'/'.join('yes' if z['hold']['target_reached'] else 'NO' for z in c['runs']) for c in sorted(result['checks'],key=lambda c:c['epsilon']))
    lower,upper=result['bulk_turning_points'];undefined=[r['epsilon'] for r in result['records'] if not r.get('three_roots')]
    text=f'''## Complete preparation scan: epsilon=30 to 40, N=50

This consolidated section supersedes the earlier drive-scan discussions. **The full vacuum-relaxation plus compensation-ramp protocol has now been tested below 35 wherever a separate intermediate root exists**, not merely its stationary stability. Existence, OBC stability, topology, and preparation success remain distinct quantities.

### Root existence and the meaning of drives below 30.16063

The uniform cubic is epsilon^2=n[(c+u n)^2+gamma^2], c=1.4998, u=-0.0004, gamma=0.5. Its turning points satisfy 3u^2n^2+4cu n+c^2+gamma^2=0, giving

**{lower['epsilon']:.8f}<epsilon<{upper['epsilon']:.8f}**

for three distinct positive bulk roots. The intermediate/high roots merge at the lower endpoint, density {lower['density']:.6f}; the intermediate/low roots merge at the upper endpoint, density {upper['density']:.6f}. Root-specific endpoint compensation makes each existing bulk root an exact uniform OBC stationary solution, including endpoints. At or beyond a root merger there is no separate intermediate target to prepare.

For epsilon=30 to the lower turning point, the intermediate-root compensation is therefore undefined. Such points are marked **no intermediate target**, not failed ramp simulations; a low-density stationary solution still exists. Cubic roots were explicitly evaluated at every requested coarse drive, including those below the existence boundary. There are {len(undefined)} sampled drives without three distinct roots. No arbitrary continuation of endpoint compensation past the root merger is introduced.

### Protocol and coverage

N=50, J=1, U=-0.0002 with U n^2 normal ordering retained, V=0, Delta=0.5, kappa=1, phi=pi/3. Drive is real and uniform in the bulk. For each drive independently, start with alpha_j=0; evolve without compensation for 1200/J. Then add fixed target-specific endpoint fields with s(t)=3(t/T)^2-2(t/T)^3, testing T=1000/J and 10000/J. Hold s=1 for another 1200/J. Endpoint additive pumps are epsilon+s J exp(+i phi) alpha0 and epsilon+s J exp(-i phi) alpha0. No root-polished state is substituted during preparation.

The combined scan covers 30-40 on a 0.05 grid, retains the prior fine central and outcome-transition samples, and adds lower-edge drives 30.1607,30.161,30.162,30.165,30.17,30.18,30.19. Coarse brackets with a change of either ramp outcome or root existence are refined to 0.01. There are {len(result['records'])} sampled drives, {len(rr)} with an intermediate target; {len(low)} eligible points below 35 have now received both full ramp tests. All uncompensated vacuum-preparation endpoints pass the existing stationary criterion: {all(r['baseline']['reached_stationary'] for r in rr)}.

### Preparation results

Spans below describe consecutive **sampled** points with the same outcome, not certified continuum windows or globally optimal preparation limits. Both the density and complex phase are checked, including the end sites.

|ramp time T (1/J)|sampled drive span|sample count|uniform intermediate target reached?|
|---:|---|---:|---|
{chr(10).join(group)}

Both ramp times succeed at every sampled eligible drive below 35, including epsilon=30.1607. No lower preparation-failure boundary is observed before the intermediate root ceases to exist. This is numerical evidence down to the nearest tested point, not proof of success arbitrarily close to or exactly at the merger. The previously observed upper preparation change between 39.61 and 39.62 persists and remains distinct from the upper root-existence boundary.

![Complete preparation scan over epsilon=30-40](../figures/intermediate_preparation_30_40.png)

![Detailed low-drive preparation scan and near-existence-boundary checks](../figures/intermediate_preparation_30_35_detail.png)

The largest uniform-target growth over eligible drives is {max(r['target']['growth'] for r in rr):.6f} J. Sampled zero-frequency winding labels are {sorted(set(r['local_winding']['winding'] for r in rr),key=str)}. These characterize the existing uniform targets, independently of which attractor the ramp reaches. Failed preparation is not evidence that the target is absent or unstable; compare final-state force, root distance, temporal variation, and BdG growth in the saved records. The analytical lower existence boundary is distinct from any finite-time preparation boundary apparent in the outcome table.

### Convergence and independent checks

Base dynamics uses DOP853, rtol=2e-10, atol=2e-11. Target success requires a raw endpoint force below 1e-3, maximum complex amplitude distance below 1e-3 from a stationary root accepted with max|f|<1e-8, and max_j|alpha_j-alpha0|<1e-3. Analytic-Jacobian root tolerance is 1e-10. Eigenmodes evolve as exp(-i omega t); positive Im(omega) means growth, with stability classification tolerance 1e-7 J.

Independent repeats use a fresh vacuum initial condition, 2400/J uncompensated preparation, and tighter rtol=2e-12, atol=2e-13; both ramp times and the 1200/J final hold are retained. Selected low drives include 30.1607,30.17,31,32,33,34 and both sides of each observed low-drive change in outcome. Existing upper-drive checks are retained. Tight outcomes (T=1000/T=10000): {checks}.

For these selected uniform targets, characteristic-polynomial/Routh root counts at 100 and 160 decimal digits check the number of growing modes. These and double-precision residuals are numerical evidence, not interval certificates. Near a saddle-node, distinct densities and required preparation precision can become small; the nearest tested drive is 30.1607, above the analytical boundary by about 6.54e-5. No ramp was attempted at the exact merger. Finite ramp tests do not establish an asymptotic adiabatic limit, global uniqueness, or the full basin geometry.

### All eligible sampled drives

|epsilon / J|intermediate target density|target OBC growth / J|zero-frequency winding|T=1000 target reached?|T=10000 target reached?|
|---:|---:|---:|---:|---|---|
{chr(10).join(detail)}

Reproduce the extension with `python src/scan_preparation_30_35.py`. Per-drive low records and tighter repeats are in data/generated/compensation_n50/preparation_30_35/; combined source paths are indexed by intermediate_preparation_30_40_index.json in the parent directory. The previous 35-40 dataset supplies unchanged upper-drive runs. All new low-drive preparations start independently from vacuum. Records retain actual ramp trajectories, final complex amplitudes, stationary profiles, force, root distances, tail variation, target/approached-state growth, and winding. No Gaussian covariance equations or fluctuation feedback enter these calculations.
'''
    p=Path('docs/uniform_topological_segment_stability.md');base=p.read_text(encoding='utf-8')
    markers=['## Complete drive scan: epsilon=35 to 40, N=50','## Complete preparation scan: epsilon=30 to 40, N=50']
    cuts=[base.index(marker) for marker in markers if marker in base];base=base[:min(cuts)].rstrip() if cuts else base.rstrip()
    p.write_text(base+'\n\n'+text,encoding='utf-8')

if __name__=='__main__':main()
