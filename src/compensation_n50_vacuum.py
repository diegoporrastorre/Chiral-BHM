"""N=50 scan, endpoint continuation, and direct vacuum reachability tests."""
from scan_homogeneous_obc import *
from test_uniform_topological_segments import precise_count

DEST=Path('data/generated/compensation_n50')
DEST.mkdir(parents=True,exist_ok=True)

def compensated(s,a0,epsilon=38.51):
    m=Chain(N=50,drive=epsilon);m.pump=m.pump.astype(complex)
    m.pump[0]+=s*np.exp(1j*m.phi)*a0;m.pump[-1]+=s*np.exp(-1j*m.phi)*a0
    return m

def trajectory(m,initial,duration=1200,rtol=2e-10,atol=2e-11):
    # Dense output samples include startup and the long-time tail.
    def rhs(t,a):
        f=m.pump+(m.detuning+m.U-.5j*m.kappa+2*m.U*abs(a)**2)*a
        f[:-1]+=m.J*np.exp(-1j*m.phi)*a[1:];f[1:]+=m.J*np.exp(1j*m.phi)*a[:-1]
        return -1j*f
    sol=solve_ivp(rhs,[0,duration],initial,method='DOP853',rtol=rtol,atol=atol,
                  t_eval=np.unique(np.r_[np.linspace(0,150,601),np.linspace(150,duration,301),np.linspace(duration-50,duration,101)]))
    assert sol.success,sol.message
    a=sol.y[:,-1];r=m.stationary(a,1e-8)
    reached=r['converged'] and max(abs(r['alpha']-a))<1e-3 and max(abs(m.force(a)))<1e-3
    rec=diagnose(m,r['alpha']) if reached else dict(density=(abs(a)**2).tolist())
    rec.update(reached_stationary=bool(reached),end_force=float(max(abs(m.force(a)))),
               root_distance=float(max(abs(r['alpha']-a))),duration=duration,rtol=rtol,atol=atol,
               tail_density_std=float(np.max(np.std(abs(sol.y[:,sol.t>=duration-50])**2,axis=1))),
               times=sol.t.tolist(),density_trajectory=(abs(sol.y)**2).T.tolist(),
               final_alpha_real=a.real.tolist(),final_alpha_imag=a.imag.tolist())
    return a,rec

def main():
    path=DEST/'scan_35_45.json'
    rows=json.loads(path.read_text()) if path.exists() else []
    for direction,grid in [('up',np.linspace(35,45,1001)),('down',np.linspace(45,35,1001))]:
        existing=[r for r in rows if r['direction']==direction]
        a=(np.array(existing[-1]['alpha_real'])+1j*np.array(existing[-1]['alpha_imag'])) if existing else np.zeros(50,complex)
        for eps in grid[len(existing):]:
            m=Chain(N=50,drive=float(eps));r=m.stationary(a,1e-8)
            if not r['converged'] or max(eigvals(m.bogoliubov(r['alpha'])).imag)>-1e-7:
                a,_,_=evolve(m,a,1200,2e-10,2e-11);r=m.stationary(a,1e-8)
            if not r['converged']:raise RuntimeError(f'Unresolved stationary solve {direction} {eps}')
            a=r['alpha'];rec=diagnose(m,a);rec.update(epsilon=round(float(eps),6),direction=direction);rows.append(rec)
            if len(rows)%200==0:path.write_text(json.dumps(rows));print('scan',direction,eps,flush=True)
        path.write_text(json.dumps(rows))
    eps=38.51;c=1.4998;u=-.0004
    densities=sorted(np.roots([u*u,2*c*u,c*c+.25,-eps*eps]).real)
    result=dict(N=50,epsilon=eps,parameters=dict(J=1,U=-.0002,V=0,delta=.5,kappa=1,phi=float(np.pi/3)),branches=[])
    for label,n in zip(['low','middle','high'],densities):
        a0=-eps/(c+u*n-.5j);a=np.full(50,a0);m=compensated(1,a0)
        assert max(abs(m.force(a)))<1e-10
        b=dict(label=label,density_target=float(n),alpha0=[float(a0.real),float(a0.imag)],
               endpoint_pumps=[[float(m.pump[j].real),float(m.pump[j].imag)] for j in [0,-1]],
               target=diagnose(m,a),precision_checks=[precise_count(50,.4998-.0008*n,.0004*n,p) for p in [100,160]])
        # Abruptly apply the fixed compensated drives from photon vacuum.
        _,b['vacuum']=trajectory(m,np.zeros(50,complex))
        _,b['vacuum_refined']=trajectory(m,np.zeros(50,complex),2400,2e-12,2e-13)
        rng=np.random.default_rng(20261104)
        _,b['perturbed_target']=trajectory(m,a+1e-4*(rng.normal(size=50)+1j*rng.normal(size=50)))
        for name in ['vacuum','vacuum_refined','perturbed_target']:
            z=b[name];af=np.array(z['final_alpha_real'])+1j*np.array(z['final_alpha_imag'])
            z['target_amplitude_error']=float(max(abs(af-a0)))
            z['target_density_error']=float(max(abs(abs(af)**2-n)))
            z['target_reached']=bool(z['reached_stationary'] and z['target_amplitude_error']<1e-3)
        # Local root continuation, no evolution fallback masquerading as continuation.
        s=1.;step=1e-7 if label=='middle' else .005;prev=None;cont=[dict(s=s,**diagnose(m,a))];fails=[]
        for iteration in range(4000):
            sn=max(0,s-step);m=compensated(sn,a0)
            initial=a if prev is None else a+(sn-s)*(a-prev[1])/(s-prev[0])
            r=m.stationary(initial,1e-8)
            if not r['converged'] or max(abs(r['alpha']-a))>1:
                fails.append(dict(s=sn,step=step,residual=r['residual']));step/=2
                if step<1e-12:break
                continue
            prev=(s,a.copy());a=r['alpha'];s=sn;cont.append(dict(s=s,**diagnose(m,a)));step=min(.005,step*1.3)
            if s==0:break
        b.update(continuation=cont,continuation_end_s=s,failed_attempts=fails)
        # Repeat the previous staircase protocol, initialized at the target root.
        a=np.full(50,a0);stairs=[]
        for s in np.linspace(1,0,51):
            m=compensated(float(s),a0);a,z=trajectory(m,a)
            z.pop('times');z.pop('density_trajectory');z['s']=float(s);stairs.append(z)
        b['staircase']=stairs
        result['branches'].append(b);(DEST/'compensation.json').write_text(json.dumps(result))
        print(label,'target growth',b['target']['growth'],'vacuum reached',b['vacuum']['target_reached'],'refined',b['vacuum_refined']['target_reached'],'continuation end',b['continuation_end_s'],flush=True)
    plot(rows,result)

def plot(rows,result):
    fig,axes=plt.subplots(2,3,figsize=(14,8),layout='constrained')
    for col,b in enumerate(result['branches']):
        z=b['vacuum'];t=np.array(z['times']);n=np.array(z['density_trajectory'])
        for j in [0,9,24,39,49]:axes[0,col].plot(t,n[:,j],label=f'j={j+1}')
        axes[0,col].axhline(b['density_target'],color='k',ls='--',label='target');axes[0,col].set(xlim=(0,150),xlabel='time Jt',ylabel='density',title=b['label']+' compensation, vacuum start');axes[0,col].legend(fontsize=7,ncol=2)
        for key,label in [('vacuum','vacuum t=1200'),('vacuum_refined','vacuum t=2400, tighter tolerance'),('perturbed_target','perturbed target')]:axes[1,col].plot(range(1,51),b[key]['density'],label=label)
        axes[1,col].axhline(b['density_target'],color='k',ls='--',label='uniform target');axes[1,col].set(xlabel='site j',ylabel='final density');axes[1,col].legend(fontsize=7)
        if b['label']!='middle':axes[1,col].set_ylim(.9*b['density_target'],1.1*b['density_target'])
        axes[1,col].ticklabel_format(axis='y',style='plain',useOffset=False)
    fig.savefig('figures/compensation_n50_vacuum.png',dpi=180)
    fig,ax=plt.subplots(2,3,figsize=(14,8),layout='constrained')
    for col,b in enumerate(result['branches']):
        rr=b['continuation'];ax[0,col].plot([z['s'] for z in rr],[z['growth'] for z in rr]);ax[0,col].axhline(0,color='k',ls='--');ax[0,col].invert_xaxis();ax[0,col].set(title=b['label']+' continuation',xlabel='compensation s',ylabel='growth / J')
        for s in [1,.98,.5,0]:
            z=min(b['staircase'],key=lambda z:abs(z['s']-s));ax[1,col].plot(range(1,51),z['density'],label=f"s={z['s']:.2f}")
        ax[1,col].legend(fontsize=8);ax[1,col].set(xlabel='site j',ylabel='density')
    fig.savefig('figures/compensation_n50_continuation.png',dpi=180)
    fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for col,rr in enumerate([rows[:1001],list(reversed(rows[1001:]))]):
        for row,key in [(0,'density'),(1,'nu')]:
            im=ax[row,col].imshow(np.array([z[key] for z in rr]).T,origin='lower',aspect='auto',extent=[34.995,45.005,.5,50.5],cmap='viridis' if row==0 else 'coolwarm');fig.colorbar(im,ax=ax[row,col],label=key);ax[row,col].set(xlabel='epsilon / J',ylabel='site j',title=['up','down'][col])
    fig.savefig('figures/compensation_n50_scan.png',dpi=180)

if __name__=='__main__':main()
