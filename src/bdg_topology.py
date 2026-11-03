"""Point-gap/SVD topology of the fixed homogeneous BdG chain.

Q=omega0 I-H(k), k increasing from -pi to pi; winding sign follows this gauge.
Uniform loss gamma=kappa/2, g is the BdG off-diagonal pairing coefficient.
"""
import numpy as np
from scipy.linalg import svd
from bulk_stability import matrix


def bloch(k,J,phi,delta,g,gamma):
    k=np.asarray(k)
    a=2*J*np.sin(phi)*np.sin(k)-1j*gamma
    d=delta+2*J*np.cos(phi)*np.cos(k)
    H=np.zeros(k.shape+(2,2),complex)
    H[...,0,0]=a+d; H[...,1,1]=a-d
    H[...,0,1]=g; H[...,1,0]=-np.conj(g)
    return H


def winding(J,phi,delta,g,gamma,omega0=0.,nk=8192):
    if abs(omega0)<1e-14 and abs(np.sin(phi))<1e-12 and abs(g)>=gamma:
        q=np.sqrt(abs(g)**2-gamma**2)
        if any(abs(target-delta)<=2*abs(J*np.cos(phi)) for target in [-q,q]):
            return dict(winding=None,raw=None,bulk_singular_gap=0.,max_phase_step=None)
    k=np.linspace(-np.pi,np.pi,nk,endpoint=False)
    Q=omega0*np.eye(2)-bloch(k,J,phi,delta,g,gamma)
    det=np.linalg.det(Q)
    increments=np.angle(np.roll(det,-1)/det)
    gap=float(np.linalg.svd(Q,compute_uv=False)[...,-1].min())
    value=float(increments.sum()/(2*np.pi))
    return dict(winding=None if gap<1e-8 or abs(increments).max()>np.pi*.9 else int(round(value)),raw=value,
                bulk_singular_gap=gap,max_phase_step=float(abs(increments).max()))


def zero_frequency_interval(J,delta,g,gamma):
    """For delta>0, 0<phi<pi/2, gamma>0; return open interval in degrees."""
    if not (J>0 and delta>0 and gamma>0): raise ValueError('Outside documented domain')
    if abs(g)<=gamma: return None
    q=np.sqrt(abs(g)**2-gamma**2)
    lo=abs(delta-q)/(2*J); hi=(delta+q)/(2*J)
    if lo>=1: return None
    return [float(np.rad2deg(np.arccos(min(1,hi)))),float(np.rad2deg(np.arccos(lo)))]


def finite_svd(N,J,phi,delta,g,gamma,omega0=0.):
    Q=omega0*np.eye(2*N)-matrix(N,J,phi,delta,g,gamma)
    u,s,vh=svd(Q)
    v=vh.conj().T[:,-1]; left=u[:,-1]
    vw=abs(v[:N])**2+abs(v[N:])**2
    uw=abs(left[:N])**2+abs(left[N:])**2
    edge=max(1,N//10)
    return dict(N=N,smin=float(s[-1]),snext=float(s[-2]),
                right_center=float(np.dot(np.arange(1,N+1),vw)),
                left_center=float(np.dot(np.arange(1,N+1),uw)),
                right_edge_weights=[float(vw[:edge].sum()),float(vw[-edge:].sum())],
                left_edge_weights=[float(uw[:edge].sum()),float(uw[-edge:].sum())],
                residual=float(np.linalg.norm(Q@v-s[-1]*left)))
