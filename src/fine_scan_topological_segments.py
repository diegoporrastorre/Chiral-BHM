"""Fine two-direction stationary scan and fixed-background segment stability."""
from scan_homogeneous_obc import *
from test_uniform_topological_segments import spectrum, precise_count

def main():
    target=OUT/'fine_scan_35_45.json'
    old=json.loads((OUT/'scan.json').read_text())
    records=json.loads(target.read_text()) if target.exists() else []
    for direction,grid,seed in [('up',np.linspace(35,45,1001),old[:161]),('down',np.linspace(45,35,1001),old[161:])]:
        existing=[r for r in records if r['direction']==direction]
        r0=existing[-1] if existing else next(r for r in seed if r['epsilon']==grid[0])
        a=np.array(r0['alpha_real'])+1j*np.array(r0['alpha_imag'])
        for eps in grid[len(existing):]:
            m=Chain(N=100,drive=float(eps));r=m.stationary(a,1e-8)
            fallback=False
            if not r['converged'] or max(eigvals(m.bogoliubov(r['alpha'])).imag)>-1e-7:
                fallback=True;a,std,_=evolve(m,a,1200,2e-10,2e-11);r=m.stationary(a,1e-8)
            if r['converged']:
                a=r['alpha'];rec=diagnose(m,a);rec['stationary']=True
            else:
                rec=dict(stationary=False,residual=r['residual'],density=(abs(a)**2).tolist())
            rec.update(epsilon=round(float(eps),6),direction=direction,evolution_fallback=fallback)
            records.append(rec)
            if len(records)%100==0:
                target.write_text(json.dumps(records));print('stationary',direction,eps,flush=True)
        target.write_text(json.dumps(records))
    # Longest contiguous candidate segment at each up-sweep point, excluding endpoints.
    tests_path=OUT/'fine_uniform_segment_tests.json'
    tests=json.loads(tests_path.read_text()) if tests_path.exists() else []
    done={r['epsilon'] for r in tests}
    for r in records[:1001]:
        eps=r['epsilon']
        if eps in done or not 35<eps<45 or not r['stationary']:continue
        ids=np.where(np.array(r['nu'])!=0)[0]
        groups=[g for g in np.split(ids,np.where(np.diff(ids)>1)[0]+1) if len(g)]
        rec=dict(epsilon=eps,tests=[])
        if groups:
            group=max(groups,key=len);density=np.array(r['density'])[group]
            rec.update(start=int(group[0]+1),end=int(group[-1]+1),length=len(group))
            for label,n in [('min',density.min()),('median',np.median(density)),('max',density.max())]:
                for length in sorted(set([len(group),100])):
                    d=.4998-.0008*n;g=.0004*n;z=spectrum(length,d,g)
                    z.update(sample=label,density=float(n),d=float(d),g=float(g));rec['tests'].append(z)
        tests.append(rec)
        if len(tests)%100==0:
            tests_path.write_text(json.dumps(tests));print('uniform tests',eps,flush=True)
    tests_path.write_text(json.dumps(tests))
    up=records[:1001];down=list(reversed(records[1001:]))
    diff=[float(max(abs(np.array(a['density'])-np.array(b['density'])))) for a,b in zip(up,down)]
    fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for col,rr in enumerate([up,down]):
        for row,key,cmap in [(0,'density','viridis'),(1,'nu','coolwarm')]:
            z=np.array([r.get(key,[np.nan]*100) for r in rr])
            im=ax[row,col].imshow(z.T,origin='lower',aspect='auto',extent=[34.995,45.005,.5,100.5],cmap=cmap)
            fig.colorbar(im,ax=ax[row,col],label='density' if row==0 else 'local winding')
            ax[row,col].set(xlabel='epsilon / J',ylabel='site j',title=['up sweep','down sweep'][col])
    fig.savefig('figures/fine_scan_35_45_maps.png',dpi=180)
    fig,ax=plt.subplots(2,1,figsize=(10,8),layout='constrained')
    for label in ['min','median','max']:
        rr=[(r,next(z for z in r['tests'] if z['N']==100 and z['sample']==label)) for r in tests if r.get('length',0)>=5]
        ax[0].plot([r['epsilon'] for r,z in rr],[z['growth'] for r,z in rr],label=label+' density')
    ax[0].axhline(0,color='k',ls='--');ax[0].legend();ax[0].set(xlabel='epsilon / J',ylabel='uniform OBC N=100 growth / J',title='Longest segment with at least five sites')
    for rr,label in [(up,'up'),(down,'down')]:ax[1].plot([r['epsilon'] for r in rr],[r.get('growth',np.nan) for r in rr],label=label)
    ax[1].axhline(0,color='k',ls='--');ax[1].legend();ax[1].set(xlabel='epsilon / J',ylabel='full inhomogeneous growth / J')
    fig.savefig('figures/fine_scan_35_45_stability.png',dpi=180)
    best=max(up,key=lambda r:sum((np.array(r.get('nu',[]))!=0)&(np.array(r.get('variation',[]))<.05)))
    fig,axes=plt.subplots(2,4,figsize=(16,7),sharex=True,sharey=True,layout='constrained')
    chosen=[min(up,key=lambda r:abs(r['epsilon']-(best['epsilon']+offset))) for offset in [-.03,-.02,-.01,0,.01,.02,.03,.04]]
    for ax,r in zip(axes.flat,chosen):
        n=np.array(r['density']);mask=np.array(r['nu'])!=0
        ax.plot(range(1,101),n);ax.plot(np.where(mask)[0]+1,n[mask],'.',color='#b05b00',ms=4)
        ax.set(title=f"epsilon={r['epsilon']:.2f}",xlabel='site j',ylabel='density',ylim=(0,6500));ax.grid(alpha=.15)
    fig.savefig('figures/fine_scan_35_45_profiles.png',dpi=180)
    flat=[(r,z) for r in tests for z in r['tests']]
    worst_r,worst=max(flat,key=lambda rz:rz[1]['growth'])
    hp=[precise_count(worst['N'],worst['d'],worst['g'],p) for p in [100,160]]
    summary=dict(points=len(records),max_residual=max(r['residual'] for r in records),
                 max_growth=max(r.get('growth',float('inf')) for r in records),failed=sum(not r['stationary'] for r in records),
                 max_direction_density_difference=max(diff),different_drives=[r['epsilon'] for r,d in zip(up,diff) if d>.01],
                 best_slow_domain_epsilon=best['epsilon'],best_slow_sites=(np.where((np.array(best['nu'])!=0)&(np.array(best['variation'])<.05))[0]+1).tolist(),
                 uniform_cases=len(flat),uniform_unstable=sum(z['growth']>1e-7 for r,z in flat),worst_uniform=dict(epsilon=worst_r['epsilon'],**worst),high_precision=hp)
    (OUT/'fine_scan_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':main()
