"""Resolve extremely small intermediate-branch boundary mismatch."""
from endpoint_compensation import *
from scipy.linalg import solve

def main():
    path=OUT/'endpoint_compensation.json';result=json.loads(path.read_text())
    for b in result['branches']:
        a0=complex(*b['alpha0']);a=np.full(100,a0);m=model(1,a0)
        forcing=np.zeros(100,complex);forcing[0]=np.exp(1j*m.phi)*a0;forcing[-1]=np.exp(-1j*m.phi)*a0
        tangent=solve(m.bogoliubov(a),-np.r_[forcing,-forcing.conj()])[:100]
        b['initial_static_susceptibility']=dict(max_dalpha_ds=float(max(abs(tangent))),
                   norm_dalpha_ds=float(np.linalg.norm(tangent)))
        if b['label']!='middle':continue
        rr=[state(m,a,1)];s=1.;step=1e-11;prev=None;fails=[]
        for iteration in range(3000):
            sn=max(0,s-step);m=model(sn,a0)
            initial=a+(sn-s)*tangent if prev is None else a+(sn-s)*(a-prev[1])/(s-prev[0])
            r=m.stationary(initial,1e-8);change=float(max(abs(r['alpha']-a)))
            if not r['converged'] or change>1:
                fails.append(dict(s=sn,step=step,residual=r['residual'],change=change));step/=2
                if step<1e-13:break
                continue
            prev=(s,a.copy());a=r['alpha'];s=sn;rr.append(state(m,a,s));step=min(.005,step*1.3)
            if s==0:break
        b['continuation']=rr;b['failed_attempts']=fails;b['continuation_end_s']=s
        print('middle refined end',s,'mismatch',1-s,'points',len(rr),'growth',rr[-1]['growth'],flush=True)
    path.write_text(json.dumps(result,indent=2));plot(result);plot_middle(result)

def plot_middle(result):
    fig,ax=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    b=result['branches'][1];rr=b['continuation']
    positive=[r for r in rr if r['s']<1]
    ax[0].plot([1-r['s'] for r in positive],[r['growth'] for r in positive],'.-');ax[0].axhline(0,color='k',ls='--')
    ax[0].set(xlabel='removed compensation 1-s',ylabel='growth / J',title='Intermediate stationary branch',xscale='log',xlim=(5e-12,1))
    ax[0].set_xticks([1e-10,1e-8,1e-6,1e-4,1e-2,1]);ax[0].minorticks_off()
    indices=[0]+[min(range(len(rr)),key=lambda i:abs((1-rr[i]['s'])-t)) for t in [1e-9,1e-7,1e-5,.01]]+[len(rr)-1]
    for index in sorted(set(indices)):
        r=rr[index];ax[1].plot(range(1,101),r['density'],label=f"1-s={1-r['s']:.3g}")
    ax[1].legend();ax[1].set(xlabel='site j',ylabel='density')
    fig.savefig('figures/endpoint_compensation_middle_refinement.png',dpi=180)

if __name__=='__main__':main()
