"""Analytic and numerical zero-frequency topology across detuning and loss, V=0."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json
import numpy as np
from bdg_topology import winding
from mean_field import Chain


def main():
    U=-.0002
    records=[]; summary=[]
    phases=np.deg2rad(np.linspace(2,178,45))  # includes 90 degrees
    for delta in [-1,0,.5,1,2,3,4]:
        for kappa in [.1,.3,1,2,4,6]:
            a=delta+U; gamma=kappa/2; local=[]
            for phi in phases:
                h=2*np.cos(phi);c=a+h
                if c<=np.sqrt(3)*gamma:continue
                rad=np.sqrt(c*c-3*gamma*gamma)
                for frac in np.linspace(.01,.99,21):
                    t=(2*c-rad)/3+frac*2*rad/3
                    n=t/(-2*U);d=a-2*t
                    eps=np.sqrt(n*((c-t)**2+gamma*gamma))
                    x0=t*t-gamma*gamma-(d+h)**2
                    xp=t*t-gamma*gamma-(d-h)**2
                    assert x0>0
                    nu=-1 if xp<0 else 0
                    records.append([delta,kappa,phi,n,eps,d,t,x0,xp,nu])
                    local.append(nu)
            summary.append(dict(delta=delta,kappa=kappa,samples=len(local),
                                topological=local.count(-1),trivial=local.count(0)))
    arr=np.array(records);checks=[]
    for nk in [2048,4096]:
        k=np.linspace(-np.pi,np.pi,nk,endpoint=False);maxstep=0.;maxerror=0.
        for start in range(0,len(arr),64):
            r=arr[start:start+64];p=r[:,2,None];d=r[:,5,None];g=r[:,6,None];gamma=r[:,1,None]/2
            det=(2*np.sin(p)*np.sin(k)-1j*gamma)**2-(d+2*np.cos(p)*np.cos(k))**2+g*g
            step=np.angle(np.roll(det,-1,axis=1)/det)
            raw=step.sum(axis=1)/(2*np.pi)
            assert np.array_equal(np.rint(raw).astype(int),r[:,9].astype(int))
            maxerror=max(maxerror,float(np.max(np.abs(raw-r[:,9]))))
            maxstep=max(maxstep,float(np.max(np.abs(step))))
        checks.append(dict(nk=nk,mismatches=0,max_integer_error=maxerror,max_phase_step=maxstep))
        print(checks[-1],flush=True)
    examples=[]
    for delta,kappa,phi,t in [(3.0002,1,np.pi/3,1.5),(.5,.1,np.pi/2,2*(.5+U)/3),(.5,1,np.pi/3,.0004*1990.520824756707)]:
        a=delta+U;c=a+2*np.cos(phi);gamma=kappa/2;n=t/(-2*U)
        eps=np.sqrt(n*((c-t)**2+gamma*gamma));alpha=-eps/(c-t-1j*gamma)
        roots=np.sort(np.roots([4*U*U,4*c*U,c*c+gamma*gamma,-eps*eps]).real)
        assert abs(roots[1]-n)<1e-6
        m=Chain(N=20,delta=delta,kappa=kappa,phi=phi,drive=eps)
        m.pump=m.pump.astype(complex);m.pump[0]+=np.exp(1j*phi)*alpha;m.pump[-1]+=np.exp(-1j*phi)*alpha
        residual=float(np.abs(m.force(np.full(20,alpha))).max());assert residual<1e-9
        tests=[winding(1,phi,a-2*t,t,gamma,nk=nk) for nk in [8192,16384]]
        assert tests[0]['winding']==tests[1]['winding']
        examples.append(dict(delta=delta,kappa=kappa,phi=float(phi),density=float(n),epsilon=float(eps),
                             roots=roots.tolist(),X0=float(t*t-gamma*gamma-(a-2*t+2*np.cos(phi))**2),
                             Xpi=float(t*t-gamma*gamma-(a-2*t-2*np.cos(phi))**2),
                             residual=residual,winding_checks=tests))
    dest=Path('data/generated/intermediate_topology_general');dest.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(dest/'scan.npz',records=arr,columns=np.array(['Delta','kappa','phi','density','epsilon','effective_detuning','pairing','X0','Xpi','winding']))
    result=dict(parameters=dict(J=1,U=U,V=0,omega_ref=0),samples=len(arr),summary=summary,checks=checks,examples=examples)
    (dest/'results.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(dict(samples=len(arr),topological=int((arr[:,9]==-1).sum()),trivial=int((arr[:,9]==0).sum()),examples=examples),indent=2))


if __name__=='__main__':main()
