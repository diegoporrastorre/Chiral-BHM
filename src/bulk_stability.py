"""Fixed-background g_c=0 BdG bulk and boundary diagnostics.

Units: hbar=a=J_unit=1. g is the off-diagonal BdG pairing, gamma=kappa/2.
No mean-field solver or covariance feedback. See derivation for limitations.
"""
import numpy as np
from numpy.polynomial import Polynomial as P
from scipy.linalg import eigvals
from scipy.optimize import linear_sum_assignment


def matrix(N, J, phi, delta, g, gamma=0., pbc=False):
    A = np.eye(N, dtype=complex)*delta
    for j in range(N if pbc else N-1):
        k = (j+1) % N
        A[k,j] = J*np.exp(1j*phi)
        A[j,k] = J*np.exp(-1j*phi)
    G = np.eye(N)*g
    return np.block([[A,G],[-G.conj(),-A.conj()]])-1j*gamma*np.eye(2*N)


def saddles(J, phi, delta, g):
    """Candidate stationary saddles; NOT a pinch/contour selection algorithm."""
    s, c = J*np.sin(phi), J*np.cos(phi)
    z = P([0,1]); u = z*z+1; w = z*z-1; q = delta*z+c*u
    poly = s*s*u*u*(q*q-abs(g)**2*z*z)+c*c*q*q*w*w
    out = []
    for b in poly.roots():
        if abs(b)<1e-10: continue
        k = -1j*np.log(complex(b))
        a = 2*s*np.sin(k); d = delta+2*c*np.cos(k)
        r = np.sqrt(complex(d*d-abs(g)**2))
        if abs(r)<1e-9: continue
        for sign in [-1,1]:
            residual = abs(2*s*np.cos(k)-sign*2*c*d*np.sin(k)/r)
            if residual < 1e-6:
                out.append((complex(k),complex(a+sign*r),float(residual)))
    return out


def local_response(times, J, phi, delta, g, gamma=0., nk=16384):
    """Infinite-chain on-site 2x2 propagator, real-k trapezoidal quadrature.

    Accurate only before cancellation of growing Fourier components dominates.
    """
    k = np.arange(nk)*2*np.pi/nk
    a = 2*J*np.sin(phi)*np.sin(k)
    d = delta+2*J*np.cos(phi)*np.cos(k)
    r = np.sqrt((d*d-abs(g)**2).astype(complex))
    vals = []
    for t in times:
        f = np.exp(-1j*a*t-gamma*t)
        co = np.cos(r*t)
        si = np.full(nk,complex(t))
        np.divide(np.sin(r*t),r,out=si,where=abs(r)>1e-12)
        block = np.array([[np.mean(f*(co-1j*d*si)),np.mean(-1j*f*g*si)],
                          [np.mean(1j*f*np.conj(g)*si),np.mean(f*(co+1j*d*si))]])
        vals.append(np.linalg.norm(block))
    return np.array(vals)


def spatial_roots(omega, J, phi, delta, g):
    """Roots z=exp(ik) of det(omega-H(k))=0, with zero uniform loss."""
    z = P([0,1]); s=J*np.sin(phi); c=J*np.cos(phi)
    q=delta*z+c*(z*z+1)
    poly=(omega*z+1j*s*(z*z-1))**2-q*q+abs(g)**2*z*z
    return poly.roots()


def pinch_candidate(k, omega, J, phi, delta, g, steps=1600, epsilon=1e-7):
    """Track the coalescing pair vertically upward in omega.

    A numerical pinch diagnostic, not a general contour-topology theorem.
    Refine steps and epsilon; a path encountering other branch points needs
    separate treatment. Returns membership in the small/large-z sets at high Im.
    """
    args=(J,phi,delta,g)
    heights=np.geomspace(epsilon,10*(abs(J)+abs(delta)+abs(g)+1),steps)
    roots=spatial_roots(omega+1j*heights[0],*args)
    pair=np.argsort(abs(roots-np.exp(1j*k)))[:2]
    worst=0.
    for height in heights[1:]:
        new=spatial_roots(omega+1j*height,*args)
        row,col=linear_sum_assignment(abs(roots[:,None]-new[None,:]))
        new=new[col]
        separation=min(abs(roots[i]-roots[j]) for i in range(4) for j in range(i))
        worst=max(worst,float(np.max(abs(new-roots))/(separation+1e-30)))
        roots=new
    membership=abs(roots[pair])<1
    return bool(membership[0]!=membership[1]),worst


def bulk_threshold(J, phi, delta, g, steps=1600, epsilon=1e-7):
    """Lossless fixed-position growth via upward pinch tracking.

    Exact endpoint expressions avoid polynomial degeneracies. Uniform loss
    shifts the returned rate by -gamma. Intended for the documented sweep.
    """
    if abs(np.sin(phi))<1e-10:
        dmin=max(0,abs(delta)-2*abs(J*np.cos(phi)))
        return float(np.sqrt(max(0,abs(g)**2-dmin**2))),[]
    if abs(np.cos(phi))<1e-10:
        return float(np.sqrt(max(0,abs(g)**2-delta**2))),[]
    candidates=[]
    for k,w,residual in saddles(J,phi,delta,g):
        if w.imag<1e-8: continue
        if any(abs(w-complex(r['omega_real'],r['omega_imag']))<1e-6 for r in candidates): continue
        pinched,step_ratio=pinch_candidate(k,w,J,phi,delta,g,steps,epsilon)
        candidates.append(dict(omega_real=w.real,omega_imag=w.imag,k_real=k.real,k_imag=k.imag,
                               pinched=pinched,step_ratio=step_ratio,residual=residual))
    return max([0.]+[r['omega_imag'] for r in candidates if r['pinched']]),candidates


def benchmark():
    """Reproducible example plus consistency checks; returns JSON-safe results."""
    from scipy.linalg import eig
    args=dict(J=1.,phi=np.pi/3,delta=.5,g=.4)
    gamma=.35
    rows=[]
    for N in [20,40,80,160,240]:
        H=matrix(N,**args)
        w,v=eig(H)
        idx=np.argmax(w.imag)
        residual=np.linalg.norm(H@v[:,idx]-w[idx]*v[:,idx])/np.linalg.norm(H)
        rows.append(dict(N=N,growth_no_loss=float(w[idx].imag),
                         growth_with_loss=float(w[idx].imag-gamma),
                         eigen_residual=float(residual)))
    N=20
    k=np.arange(N)*2*np.pi/N
    d=args['delta']+2*np.cos(args['phi'])*np.cos(k)
    a=2*np.sin(args['phi'])*np.sin(k)
    r=np.sqrt((d*d-args['g']**2).astype(complex))
    expected=np.r_[a+r,a-r]
    actual=eigvals(matrix(N,**args,pbc=True))
    pbc_error=max(min(abs(actual-x)) for x in expected)
    assert pbc_error<1e-10
    shifted=eigvals(matrix(N,**args,gamma=gamma))
    original=eigvals(matrix(N,**args))
    shift_error=max(min(abs(shifted-(x-1j*gamma))) for x in original)
    assert shift_error<1e-10
    # Exactly solvable phi=pi/2: spatial hopping commutes with onsite Nambu block.
    special=eigvals(matrix(N,1,np.pi/2,.2,.4))
    special_error=abs(special.imag.max()-np.sqrt(.4**2-.2**2))
    assert special_error<1e-10
    times=np.linspace(20,60,801)
    fits=[]
    ref=None
    for nk in [4096,8192,16384]:
        y=local_response(times,**args,nk=nk)
        fits.append(dict(nk=nk,fit_20_60=float(np.polyfit(times,np.log(y*np.sqrt(times)),1)[0]),
                         fit_40_60=float(np.polyfit(times[times>=40],np.log(y[times>=40]*np.sqrt(times[times>=40])),1)[0]),
                         relative_difference=None if ref is None else float(np.max(abs(y-ref)/y))))
        ref=y
    return dict(parameters={**args,'gamma':gamma,'kappa':2*gamma},obc=rows,
                pbc_growth=.4-gamma,local_response_fits=fits,
                candidates=[dict(k_real=k.real,k_imag=k.imag,omega_real=w.real,omega_imag=w.imag,residual=e)
                            for k,w,e in saddles(**args)],
                checks=dict(pbc_error=float(pbc_error),loss_shift_error=float(shift_error),special_case_error=float(special_error)))


if __name__ == '__main__':
    for delta, g, phi in [(0,.4,np.pi/3),(.5,.4,np.pi/3),(0,.4,np.pi/2),(0,1.2,np.pi/3)]:
        print('PARAMETERS',delta,g,phi,flush=True)
        print('SADDLES',[(round(w.real,6),round(w.imag,6)) for k,w,e in saddles(1,phi,delta,g)],flush=True)
        for N in [20,40,80,160]:
            print('OBC',N,float(eigvals(matrix(N,1,phi,delta,g)).imag.max()),flush=True)
        ts=np.linspace(10,30,201)
        response=local_response(ts,1,phi,delta,g)
        print('LOCAL',[(float(t),float(v)) for t,v in zip(ts[::50],response[::50])],flush=True)
