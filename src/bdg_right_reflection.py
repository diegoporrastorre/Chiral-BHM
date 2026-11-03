"""Right-layer modal reflection, measured at the last uniform-lead site.

Amplitudes use unit Euclidean norm of the two bulk mode vectors, not flux.
In unstable systems this is a formal frequency-domain scattering diagnostic.
"""
from decimal import localcontext
import numpy as np
from scipy.linalg import solve
from bdg_right_absorber import generator


def reflection_decimal(length,family,eta,kappa,omegas,precision=100,transparent_reference=False,n=100):
    from decimal_complex_banded import C,solve_banded
    out=[]
    with localcontext() as ctx:
        ctx.prec=precision
        a=generator(length+1,length,family,eta)
        for omega in omegas:
            u=C(str(kappa))/2-C(0,str(omega))
            disc=(u*u+C('2.61')).sqrt()*(u*u+C('3.41')).sqrt()
            zr=(u*u-C('4.99')-disc)/(4*(u+1))
            zl=-(u-1)/((u+1)*zr)
            pr=C('.1')/(u-1-2*zr); qr=C(1)
            nr=(pr.norm2()+qr.norm2()).sqrt(); pr=pr/nr; qr=qr/nr
            ql=-C('.1')*zl/((u+1)*zl+2)
            pl=C(1)/(1+ql.norm2()).sqrt()
            chi=-C('.1')/((u+1)*zl+2)
            m=[[-C(v) for v in row] for row in a]
            for i in range(len(a)): m[i][i]=m[i][i]+u
            if transparent_reference:
                assert family=='loss' and eta==0
                m[-2][-2]=m[-2][-2]-2*zr
            m[1][0]=m[1][0]+2*chi
            rhs=[C(0) for _ in a]; rhs[1]=-2*(qr/zr-pr*chi)
            sol=solve_banded(m,rhs)
            rr=(sol[0]-pr)/pl
            residual=max((sum((m[i][j]*sol[j] for j in range(max(0,i-2),min(len(a),i+3))),C(0))-rhs[i]).magnitude() for i in range(len(a)))
            # A_M/B_M needed to satisfy the sharp LEFT boundary q_0=0.
            # Normalizations here are the same as in the reflection amplitude.
            qlnorm=ql*pl
            ratio=zl/zr; power=C(1)
            for _ in range(n-length): power=power*ratio
            critical=-qr/qlnorm*power
            out.append(dict(omega=float(omega),abs_r=float(rr.magnitude()),r_real=float(rr.re),r_imag=float(rr.im),
                            unit_loop_reflection_scale=float(critical.magnitude()),
                            absolute_residual=float(residual),precision=precision))
    return out


def reflection(length,family,eta,kappa,omegas):
    with localcontext() as ctx:
        ctx.prec=40
        a=np.array([[float(v) for v in row] for row in generator(length+1,length,family,eta)])
    out=[]
    for omega in omegas:
        u=kappa/2-1j*omega
        disc=np.sqrt(u*u+2.61)*np.sqrt(u*u+3.41)
        aa=u*u-4.99
        zr=(aa-disc)/(4*(u+1))
        zl=-(u-1)/((u+1)*zr)
        vr=np.array([.1/(u-1-2*zr),1],complex); vr/=np.linalg.norm(vr)
        vl=np.array([1,-.1*zl/((u+1)*zl+2)],complex); vl/=np.linalg.norm(vl)
        chi=-.1/((u+1)*zl+2)
        m=u*np.eye(len(a),dtype=complex)-a
        m[1,0]+=2*chi
        rhs=np.zeros(len(a),complex)
        rhs[1]=-2*(vr[1]/zr-vr[0]*chi)
        sol=solve(m,rhs)
        rr=(sol[0]-vr[0])/vl[0]
        residual=np.linalg.norm(m@sol-rhs)/(np.linalg.norm(m)*np.linalg.norm(sol)+np.linalg.norm(rhs))
        out.append(dict(omega=float(omega),abs_r=float(abs(rr)),r_real=float(rr.real),r_imag=float(rr.imag),
                        solve_residual=float(residual)))
    return out


if __name__=='__main__':
    for l,f,e in [(0,'quintic',0),(40,'quintic',0),(80,'quintic',0),(90,'linear',0),(80,'loss',3)]:
        r=reflection(l,f,e,2.1,[0,.5,1,2])
        print(l,f,e,[(x['omega'],x['abs_r']) for x in r])
