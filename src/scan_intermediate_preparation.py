"""Drive dependence of vacuum-relaxation plus smooth compensation protocol."""
from adiabatic_compensation_n50 import *
from bdg_topology import winding

def compact(run):
    # Preserve ramp samples and final hold diagnostics, omit repeated long hold history.
    run['hold'].pop('times',None);run['hold'].pop('density_trajectory',None)
    return run

def main():
    path=DEST/'intermediate_preparation_scan.json'
    result=json.loads(path.read_text()) if path.exists() else dict(N=50,records=[],checks=[])
    grid=sorted(set(np.round(np.r_[np.arange(38.40,38.601,.01),[38.25,38.30,38.35,38.65,38.70,38.75]],6)))
    done={r['epsilon'] for r in result['records']}
    for eps in grid:
        if eps in done:continue
        c=1.4998;u=-.0004
        roots=np.roots([u*u,2*c*u,c*c+.25,-eps**2])
        densities=sorted(x.real for x in roots if abs(x.imag)<1e-7 and x.real>=0)
        if len(densities)!=3:
            result['records'].append(dict(epsilon=float(eps),three_roots=False));continue
        n=densities[1];a0=-eps/(c+u*n-.5j)
        m=compensated(1,a0,float(eps));target=diagnose(m,np.full(50,a0))
        a,base=trajectory(Chain(N=50,drive=float(eps)),np.zeros(50,complex),1200)
        rec=dict(epsilon=float(eps),three_roots=True,density_target=float(n),alpha0=[float(a0.real),float(a0.imag)],
                 target=target,local_winding=winding(1,np.pi/3,.4998-.0008*n,.0004*n,.5),baseline=base,runs=[])
        rec['baseline'].pop('times');rec['baseline'].pop('density_trajectory')
        for T in [1000,10000]:
            z=compact(ramp(a.copy(),a0,T,epsilon=float(eps)))
            z['hold']['max_relative_density_error']=float(max(abs(np.array(z['hold']['density'])-n))/n)
            rec['runs'].append(z)
        result['records'].append(rec);path.write_text(json.dumps(result))
        print(eps,n,target['growth'],[z['hold']['target_reached'] for z in rec['runs']],flush=True)
    if not result['checks']:
        # Independent tighter complete protocols at endpoints, centre, and all slower-ramp failures.
        selected={38.25,38.51,38.75}
        selected.update(r['epsilon'] for r in result['records'] if r.get('three_roots') and not r['runs'][-1]['hold']['target_reached'])
        for eps in sorted(selected):
            rec=next(r for r in result['records'] if r['epsilon']==eps)
            if not rec.get('three_roots'):continue
            a0=complex(*rec['alpha0'])
            a,base=trajectory(Chain(N=50,drive=eps),np.zeros(50,complex),2400,2e-12,2e-13)
            z=compact(ramp(a,a0,10000,2e-12,2e-13,epsilon=eps))
            check=dict(epsilon=eps,run=z,baseline_end_force=base['end_force'],
                       precision_checks=[precise_count(50,.4998-.0008*rec['density_target'],.0004*rec['density_target'],p) for p in [100,160]])
            result['checks'].append(check);path.write_text(json.dumps(result));print('tight',eps,z['hold']['target_reached'],flush=True)
    plot_scan(result);report_scan(result)

def plot_scan(result):
    rr=sorted([r for r in result['records'] if r.get('three_roots')],key=lambda r:r['epsilon'])
    fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for i,T in enumerate([1000,10000]):
        error=[r['runs'][i]['hold']['target_amplitude_error'] for r in rr]
        ax[0,0].semilogy([r['epsilon'] for r in rr],error,'-o',ms=4,label=f'T={T}')
        ax[0,1].plot([r['epsilon'] for r in rr],[np.mean(r['runs'][i]['hold']['density']) for r in rr],'-o',ms=4,label=f'T={T}')
    ax[0,0].axhline(1e-3,color='k',ls='--',label='target threshold');ax[0,0].legend();ax[0,0].set(xlabel='epsilon / J',ylabel='max complex amplitude error')
    ax[0,1].plot([r['epsilon'] for r in rr],[r['density_target'] for r in rr],'k--',label='uniform intermediate target');ax[0,1].legend();ax[0,1].set(xlabel='epsilon / J',ylabel='final mean density')
    ax[1,0].plot([r['epsilon'] for r in rr],[r['target']['growth'] for r in rr],'-o');ax[1,0].axhline(0,color='k',ls='--');ax[1,0].set(xlabel='epsilon / J',ylabel='uniform target OBC growth / J')
    for eps in [38.25,38.40,38.51,38.60,38.75]:
        r=next(r for r in rr if r['epsilon']==eps);ax[1,1].plot(range(1,51),r['runs'][-1]['hold']['density'],label=f'epsilon={eps:.2f}')
    ax[1,1].legend();ax[1,1].set(xlabel='site j',ylabel='final density, T=10000')
    fig.savefig('figures/intermediate_preparation_epsilon_scan.png',dpi=180)

def report_scan(result):
    rr=sorted([r for r in result['records'] if r.get('three_roots')],key=lambda r:r['epsilon'])
    lines=[]
    for r in rr:
        lines.append(f"|{r['epsilon']:.2f}|{r['density_target']:.2f}|{r['target']['growth']:.5f}|{r['local_winding']['winding']}|"+'|'.join('yes' if z['hold']['target_reached'] else 'NO' for z in r['runs'])+'|')
    successes=[[r['epsilon'] for r in rr if r['runs'][i]['hold']['target_reached']] for i in range(2)]
    text=f'''## Drive dependence of intermediate-state preparation at N=50

We repeat the complete **vacuum -> uncompensated stationary state -> smooth intermediate-root compensation -> final hold** protocol at {len(rr)} drives. The central window 38.40-38.60 is sampled every 0.01; shoulder points are 38.25,38.30,38.35,38.65,38.70,38.75. Each drive is independent: its preparation starts from vacuum, not from the previous drive's result. We recompute all three bulk roots, choose the intermediate density, and calculate its complex endpoint compensation at that drive. N=50 and all other model parameters remain unchanged.

The uncompensated preparation lasts 1200/J. Use s(t)=3(t/T)^2-2(t/T)^3 during the compensation ramp, with T=1000/J or 10000/J. Hold the final compensated Hamiltonian for another 1200/J. This scan tests finite smooth ramps, not a mathematically established adiabatic limit.

|epsilon / J|intermediate target density|target OBC growth / J|zero-frequency winding|reached with T=1000/J?|reached with T=10000/J?|
|---:|---:|---:|---:|---|---|
{chr(10).join(lines)}

The target is reached at {len(successes[0])}/{len(rr)} drives for T=1000/J and {len(successes[1])}/{len(rr)} drives for T=10000/J. These are sampled outcomes, not continuous success intervals or proven thresholds. A nonzero local winding characterizes the uniform target; stability and dynamical preparation are assessed separately.

![Drive-dependent preparation success, target stability, and final profiles](../figures/intermediate_preparation_epsilon_scan.png)

Success requires the final raw integration endpoint to lie within 1e-3 in maximum complex amplitude of the intended uniform alpha0, including all edge sites, and to pass the existing stationary convergence checks. Thus density agreement or convergence of a separate root solve alone cannot establish success. The final holds retain residuals, temporal variation, complex amplitudes, target error, density profiles, and BdG growth for the approached state.

For independent tolerance checks we rerun the full protocol at both shoulder endpoints, epsilon=38.51, and every failed T=10000/J case, with a 2400/J initial relaxation and DOP853 rtol=2e-12, atol=2e-13 instead of 2e-10, 2e-11. The final ramp and hold durations remain 10000/J and 1200/J. Check outcomes: {', '.join(str(c['epsilon'])+': '+('target reached' if c['run']['hold']['target_reached'] else 'target not reached') for c in result['checks'])}. For these target operators, 100/160-digit Routh counts are saved alongside the repeat, independently checking the number of growing eigenmodes. Such checks do not certify untested continuum parameter values or global attraction.

Reproduce with `python src/scan_intermediate_preparation.py`. The scan checkpoints per drive and uses the existing reusable mean-field solver and trajectory integrator. Data are in data/generated/compensation_n50/intermediate_preparation_scan.json, including target/background parameters, endpoint amplitudes through alpha0, baseline convergence, sampled ramp profiles, final hold diagnostics, target winding, and tighter precision checks. No fluctuation corrections feed into the mean field. Interpret failed preparation as a result of this protocol, not absence of the uniform stationary state.
'''
    p=Path('docs/uniform_topological_segment_stability.md');base=p.read_text(encoding='utf-8').split('## Drive dependence of intermediate-state preparation at N=50')[0].rstrip();p.write_text(base+'\n\n'+text,encoding='utf-8')

if __name__=='__main__':main()
