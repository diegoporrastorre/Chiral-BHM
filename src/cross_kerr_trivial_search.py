"""Search compensated cross-Kerr intermediate states for stable trivial examples."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eig, eigvals
from mean_field import Chain

U=-.0002

def symbol(k,n,phi,delta,kappa,V):
    drift=2*np.sin(phi)*np.sin(k)-.5j*kappa
    d=delta+U+(4*U+2*V)*n+2*(np.cos(phi)+V*n)*np.cos(k)
    b=2*n*(U+V*np.cos(k))
    return drift,d,b

def operator(N,n,phi,delta,kappa,V):
    C=np.eye(N,k=1)+np.eye(N,k=-1)
    T=np.exp(-1j*phi)*np.eye(N,k=1)+np.exp(1j*phi)*np.eye(N,k=-1)
    A=T+V*n*C+np.diag(np.full(N,delta+U+4*U*n-.5j*kappa)+V*n*C.sum(axis=1))
    B=2*U*n*np.eye(N)+V*n*C
    return np.block([[A,B],[-B.conj(),-A.conj()]])

def check(r,N):
    n,p,delta,kappa,V=[r[x] for x in ['density','phi','delta','kappa','V']]
    L=operator(N,n,p,delta,kappa,V)
    w,vl,vr=eig(L,left=True,right=True);ix=np.argmax(w.imag)
    condition=1/abs(np.vdot(vl[:,ix],vr[:,ix]))
    eps=r['epsilon'];alpha=-eps/(delta+U+2*np.cos(p)+2*(U+V)*n-.5j*kappa)
    model=Chain(N=N,drive=eps,delta=delta,kappa=kappa,U=U,V=V,phi=p)
    model.pump=model.pump.astype(complex)
    model.pump[0]+=(np.exp(1j*p)+V*n)*alpha
    model.pump[-1]+=(np.exp(-1j*p)+V*n)*alpha
    residual=float(np.abs(model.force(np.full(N,alpha))).max());assert residual<1e-8
    actual=model.bogoliubov(np.full(N,alpha));wa=eigvals(actual)
    assert abs(wa.imag.max()-w.imag.max())<1e-7
    numerical=[]
    for nk in [8192,16384]:
        k=np.linspace(-np.pi,np.pi,nk,endpoint=False);a,d,b=symbol(k,n,p,delta,kappa,V)
        H=np.zeros((nk,2,2),complex);H[:,0,0]=a+d;H[:,1,1]=a-d;H[:,0,1]=b;H[:,1,0]=-b
        det=np.linalg.det(-H);steps=np.angle(np.roll(det,-1)/det)
        raw=float(steps.sum()/(2*np.pi));assert abs(raw)<1e-8
        numerical.append(dict(nk=nk,winding=raw,gap=float(np.linalg.svd(H,compute_uv=False)[:,-1].min()),maxstep=float(abs(steps).max())))
    return dict(N=N,growth=float(w[ix].imag),condition=float(condition),residual=residual,
                eigen_residual=float(np.linalg.norm(L@vr[:,ix]-w[ix]*vr[:,ix])),winding=numerical)

def main():
    rows=[];candidates=[];total=0;bulk=[]
    for ratio in [1,.5,2,4,-.5]:
        V=ratio*U
        for delta in [0,.5,1,3]:
            for kappa in [.1,1,2]:
                for p in np.deg2rad(np.linspace(5,175,18)):
                    c=delta+U+2*np.cos(p);gamma=kappa/2
                    if c<=np.sqrt(3)*gamma:continue
                    rad=np.sqrt(c*c-3*gamma*gamma)
                    for frac in [1e-5,.001,.01,.1,.3,.5,.7,.9,.99,.999,1-1e-5]:
                        t=(2*c-rad)/3+2*rad*frac/3;n=t/(-2*(U+V))
                        eps=np.sqrt(n*((c-t)**2+gamma*gamma));total+=1
                        av,d,b=symbol(np.array([0,np.pi]),n,p,delta,kappa,V)
                        x=(av*av-d*d+b*b).real
                        assert x[0]>0
                        bulk.append([ratio,delta,kappa,p,n,x[0],x[1]])
                        if x[1]<=1e-8:continue
                        growth=float(eigvals(operator(20,n,p,delta,kappa,V),check_finite=False).imag.max())
                        r=dict(V=V,ratio=ratio,delta=delta,kappa=kappa,phi=float(p),density=float(n),epsilon=float(eps),X0=float(x[0]),Xpi=float(x[1]),growth20=growth)
                        r['growth50']=float(eigvals(operator(50,n,p,delta,kappa,V),check_finite=False).imag.max())
                        rows.append(r)
                        if min(growth,r['growth50']) < -1e-6:
                            r['check50']=check(r,50);candidates.append(r)
        print('ratio',ratio,'trivial so far',len(rows),'stable20',len(candidates),flush=True)
    out=Path('data/generated/cross_kerr_trivial');out.mkdir(parents=True,exist_ok=True)
    numerical=[];bb=np.array(bulk)
    for nk in [2048,4096]:
        k=np.linspace(-np.pi,np.pi,nk,endpoint=False);worst=0.
        for start in range(0,len(bb),64):
            r=bb[start:start+64];ratio,delta,kappa,p,n=[r[:,j,None] for j in range(5)]
            av,d,b=symbol(k,n,p,delta,kappa,ratio*U)
            det=av*av-d*d+b*b;step=np.angle(np.roll(det,-1,axis=1)/det)
            raw=step.sum(axis=1)/(2*np.pi);expected=np.where(r[:,6]<0,-1,0)
            assert np.array_equal(np.rint(raw),expected)
            worst=max(worst,float(abs(step).max()))
        numerical.append(dict(nk=nk,mismatches=0,max_phase_step=worst))
    result=dict(total=total,trivial=len(rows),stable20=sum(r['growth20'] < -1e-7 for r in rows),stable50=sum(r['growth50'] < -1e-7 for r in rows),rows=rows,candidates=candidates,winding_checks=numerical,
                by_ratio=[dict(ratio=q,total=int((bb[:,0]==q).sum()),trivial=sum(r['ratio']==q for r in rows),min_growth50=min([r['growth50'] for r in rows if r['ratio']==q],default=None)) for q in [1,.5,2,4,-.5]])
    np.savez_compressed(out/'bulk.npz',records=bb)
    (out/'search.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:result[k] for k in ['total','trivial','stable20','stable50']},indent=2))
    print(json.dumps(sorted(candidates,key=lambda r:r['check50']['growth'])[:3],indent=2))

if __name__=='__main__':main()
