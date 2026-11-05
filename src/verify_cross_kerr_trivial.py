"""High-precision independent growing-mode counts for cross-Kerr counterexamples."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json
from pathlib import Path
from decimal import Decimal as D,localcontext
import numpy as np
from scipy.linalg import eigvals
from banded_stability_decimal import characteristic_banded,count_right
from cross_kerr_trivial_search import operator,check

def decimal_generator(r,precision,N=50):
    """Compensated OBC time generator; retain degree-one endpoint shifts.

    Phase is specified in degrees to ten decimal places, as in the original
    verification. The density is recomputed from the intermediate fraction.
    """
    with localcontext() as ctx:
        ctx.prec=precision+20;tol=D(10)**(-precision-10)
        def atan(x):
            term=total=x;k=1
            while abs(term)>tol:
                term=-term*x*x*(2*k-1)/(2*k+1);total+=term;k+=1
            return total
        pi=16*atan(D(1)/5)-4*atan(D(1)/239)
        ph=pi*D(str(round(np.rad2deg(r['phi']),10)))/180
        cs=term=D(1);k=1
        while abs(term)>tol:term=-term*ph*ph/((2*k-1)*(2*k));cs+=term;k+=1
        sn=term=ph;k=1
        while abs(term)>tol:term=-term*ph*ph/((2*k)*(2*k+1));sn+=term;k+=1
        u=D('-.0002');v=u*D(str(r['ratio']));gamma=D(str(r['kappa']))/2;a=D(str(r['delta']))+u
        c=a+2*cs;rad=(c*c-3*gamma*gamma).sqrt()
        frac=D(str(r['fraction']));t=(2*c-rad)/3+2*rad*frac/3;n=t/(-2*(u+v))
        ctx.prec=precision;gs=-2*u*n;gc=-v*n;ar=cs+v*n
        A=[[D(0)]*(2*N) for _ in range(2*N)]
        for j in range(N):
            z=1 if j in [0,N-1] else 2;dd=a+4*u*n+z*v*n
            A[2*j][2*j]=gs-gamma;A[2*j+1][2*j+1]=-gs-gamma
            A[2*j][2*j+1]=dd;A[2*j+1][2*j]=-dd
            if j<N-1:
                for p,q,im in [(j,j+1,-sn),(j+1,j,sn)]:
                    A[2*p][2*q]=im+gc;A[2*p+1][2*q+1]=im-gc
                    A[2*p][2*q+1]=ar;A[2*p+1][2*q]=-ar
        return A,n,ph

def precise(r,precision,N=50):
    with localcontext() as ctx:
        ctx.prec=precision
        A,n,ph=decimal_generator(r,precision,N)
        v=D('-.0002')*D(str(r['ratio']))
        coeff=characteristic_banded(A,3)
        counts={str(shift):count_right(coeff,D(str(shift))) for shift in [-1e-4,0,1e-4]}
        assert all(vv==0 for vv in counts.values())
        # Check independent quadrature construction against the Nambu spectrum.
        spectrum=eigvals(np.array(A,dtype=complex))
        direct=eigvals(operator(N,float(n),float(ph),r['delta'],r['kappa'],float(v)))
        assert abs(spectrum.real.max()-direct.imag.max())<1e-8
        return dict(precision=precision,density=str(n),counts=counts)

def main():
    dest=Path('data/generated/cross_kerr_trivial')
    result=json.loads((dest/'search.json').read_text())
    chosen=[]
    for r in result['candidates']:
        if r['growth50']>=-1e-7:continue
        c=r['delta']-.0002+2*np.cos(r['phi']);gam=r['kappa']/2
        rad=np.sqrt(c*c-3*gam*gam);t=-2*(-.0002+r['V'])*r['density']
        frac=(t-(2*c-rad)/3)/(2*rad/3)
        r['fraction']=min([1e-5,.001,.01,.1,.3,.5,.7,.9,.99,.999,.99999],key=lambda f:abs(f-frac))
        r['precision_checks']=[precise(r,p) for p in [100,160]]
        r['size_checks']=[check(r,N) for N in [8,20,50,80]]
        chosen.append(r);print('verified',r['ratio'],r['fraction'],r['growth50'],flush=True)
    (dest/'verified.json').write_text(json.dumps(chosen,indent=2))
    print(json.dumps([dict(ratio=r['ratio'],fraction=r['fraction'],sizes=[(s['N'],s['growth']) for s in r['size_checks']]) for r in chosen],indent=2))

if __name__=='__main__':main()
