"""Remove fixed endpoint drives at epsilon=38.51: roots versus dynamics."""
from scan_homogeneous_obc import *
from scipy.linalg import eig

def model(s,a0):
    m=Chain(N=100,drive=38.51)
    m.pump=m.pump.astype(complex)
    m.pump[0]+=s*np.exp(1j*m.phi)*a0
    m.pump[-1]+=s*np.exp(-1j*m.phi)*a0
    return m

def state(m,a,s):
    r=diagnose(m,a);r['s']=float(s)
    return r

def main():
    eps=38.51;c=1.4998;u=-.0004
    densities=sorted(np.roots([u*u,2*c*u,c*c+.25,-eps*eps]).real)
    result=dict(epsilon=eps,N=100,parameters=dict(J=1,U=-.0002,V=0,delta=.5,kappa=1,phi=float(np.pi/3)),branches=[])
    target=OUT/'endpoint_compensation.json'
    for label,n in zip(['low','middle','high'],densities):
        a0=-eps/(c+u*n-.5j);a=np.full(100,a0)
        m=model(1,a0)
        assert max(abs(m.force(a)))<1e-10
        assert np.max(abs(m.bogoliubov(a)-Chain(N=100,drive=eps).bogoliubov(a)))==0
        branch=dict(label=label,uniform_density=float(n),alpha0=[float(a0.real),float(a0.imag)],
                    left_compensation=[float((np.exp(1j*m.phi)*a0).real),float((np.exp(1j*m.phi)*a0).imag)],
                    right_compensation=[float((np.exp(-1j*m.phi)*a0).real),float((np.exp(-1j*m.phi)*a0).imag)],
                    continuation=[state(m,a,1)],failed_attempts=[])
        s=1.;step=.005;previous=None
        while s>1e-12:
            sn=max(0,s-step);m=model(sn,a0)
            initial=a if previous is None else a+(sn-s)*(a-previous[1])/(s-previous[0])
            r=m.stationary(initial,1e-8)
            # Reject large Newton jumps; this is local continuation, not branch discovery.
            change=float(max(abs(r['alpha']-a)))
            if not r['converged'] or change>2:
                branch['failed_attempts'].append(dict(s=sn,residual=r['residual'],change=change,step=step))
                step/=2
                if step<1e-7:break
                continue
            previous=(s,a.copy());a=r['alpha'];s=sn
            branch['continuation'].append(state(m,a,s))
            step=min(.005,step*1.3)
        branch['continuation_end_s']=s
        print(label,'continuation end',s,'points',len(branch['continuation']),flush=True)
        # Actual staircase drive protocol, kept distinct from Newton continuation.
        a=np.full(100,a0);dynamics=[]
        for s in np.linspace(1,0,51):
            m=model(float(s),a0)
            a,std,sol=evolve(m,a,duration=1200,rtol=2e-10,atol=2e-11)
            endforce=float(max(abs(m.force(a))))
            r=m.stationary(a,1e-8)
            # A root counts as reached only if dynamics finishes close to it.
            distance=float(max(abs(r['alpha']-a)))
            reached=r['converged'] and distance<1e-3 and endforce<1e-3
            rec=state(m,r['alpha'],s) if reached else dict(s=float(s),density=(abs(a)**2).tolist())
            rec.update(reached_stationary=bool(reached),temporal_density_std=std,end_force=endforce,root_distance=distance)
            dynamics.append(rec)
        branch['dynamics']=dynamics
        result['branches'].append(branch);target.write_text(json.dumps(result,indent=2))
        print(label,'dynamics complete',flush=True)
    reference=next(r for r in json.loads((OUT/'fine_scan_35_45.json').read_text()) if r['epsilon']==38.51 and r['direction']=='up')
    for b in result['branches']:
        b['final_density_difference_from_original']=float(max(abs(np.array(b['dynamics'][-1]['density'])-np.array(reference['density']))))
    target.write_text(json.dumps(result,indent=2))
    plot(result)

def plot(result):
    fig,ax=plt.subplots(2,3,figsize=(14,8),layout='constrained')
    for col,b in enumerate(result['branches']):
        rr=b['continuation'];ss=[r['s'] for r in rr]
        ax[0,col].plot(ss,[r['growth'] for r in rr],label='stationary continuation')
        dd=[r for r in b['dynamics'] if r['reached_stationary']]
        ax[0,col].plot([r['s'] for r in dd],[r['growth'] for r in dd],'.',label='staircase dynamics')
        ax[0,col].axhline(0,color='k',ls='--');ax[0,col].invert_xaxis();ax[0,col].legend(fontsize=8)
        ax[0,col].set(title=b['label']+' bulk root',xlabel='compensation s',ylabel='growth / J')
        for s in [1,.98,.9,.5,0]:
            r=min(b['dynamics'],key=lambda r:abs(r['s']-s))
            ax[1,col].plot(range(1,101),r['density'],label=f"s={r['s']:.2f}")
        ax[1,col].legend(fontsize=8);ax[1,col].set(xlabel='site j',ylabel='density')
    fig.savefig('figures/endpoint_compensation_paths.png',dpi=180)
    fig,ax=plt.subplots(1,3,figsize=(14,4),layout='constrained')
    for axes,b in zip(ax,result['branches']):
        z=np.array([r['density'] for r in b['dynamics']])
        im=axes.imshow(z.T,origin='lower',aspect='auto',extent=[1.01,-.01,.5,100.5],vmin=0,vmax=6500)
        axes.set(title=b['label']+' initial bulk root',xlabel='compensation s',ylabel='site j');fig.colorbar(im,ax=axes,label='density')
    fig.savefig('figures/endpoint_compensation_density_maps.png',dpi=180)

if __name__=='__main__':main()
