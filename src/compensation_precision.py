"""Decimal banded stationary derivatives and nonlinear endpoint-error checks."""
from decimal import Decimal,localcontext
import numpy as np

def solve_real_banded(a,rhs,width=3):
    """Partial-pivot banded elimination with several real Decimal RHS columns."""
    a=[list(row) for row in a];b=[list(row) for row in rhs];n=len(a);nr=len(b[0])
    for k in range(n):
        pivot=max(range(k,min(n,k+width+1)),key=lambda i:abs(a[i][k]))
        if not a[pivot][k]:raise ArithmeticError('Singular stationary Jacobian')
        if pivot!=k:a[k],a[pivot]=a[pivot],a[k];b[k],b[pivot]=b[pivot],b[k]
        for i in range(k+1,min(n,k+width+1)):
            factor=a[i][k]/a[k][k]
            for j in range(k+1,min(n,k+2*width+1)):a[i][j]-=factor*a[k][j]
            for j in range(nr):b[i][j]-=factor*b[k][j]
            a[i][k]=Decimal(0)
    x=[[Decimal(0)]*nr for _ in range(n)]
    for i in range(n-1,-1,-1):
        for col in range(nr):
            value=b[i][col]
            for j in range(i+1,min(n,i+2*width+1)):value-=a[i][j]*x[j][col]
            x[i][col]=value/a[i][i]
    return x

def constants():
    D=Decimal;n=D('1990.52082475671');c=D('1.4998');beta=D('-.0004');e=D('38.51');gamma=D('.5')
    for _ in range(30):
        step=(n*((c+beta*n)**2+gamma*gamma)-e*e)/((c+beta*n)*(c+3*beta*n)+gamma*gamma);n-=step
        if abs(step)<D(10)**(-60):break
    h=c+beta*n;den=h*h+gamma*gamma;ar=-e*h/den;ai=-e*gamma/den
    co=D('.5');si=D(3).sqrt()/2
    ql=(co*ar-si*ai,si*ar+co*ai);qr=(co*ar+si*ai,-si*ar+co*ai)
    return n,ar,ai,co,si,ql,qr

def jacobian(x,y,co,si):
    D=Decimal;N=len(x);a=[[D(0)]*(2*N) for _ in range(2*N)];U=D('-.0002');det=D('.4998');gamma=D('.5')
    for j in range(N):
        a[2*j][2*j]=det+2*U*(3*x[j]*x[j]+y[j]*y[j])
        a[2*j+1][2*j+1]=det+2*U*(x[j]*x[j]+3*y[j]*y[j])
        a[2*j][2*j+1]=gamma+4*U*x[j]*y[j]
        a[2*j+1][2*j]=-gamma+4*U*x[j]*y[j]
        if j<N-1:
            for p,q,im in [(j,j+1,-si),(j+1,j,si)]:
                a[2*p][2*q]=co;a[2*p][2*q+1]=-im
                a[2*p+1][2*q]=im;a[2*p+1][2*q+1]=co
    return a

def derivatives(N,precision=100):
    with localcontext() as ctx:
        ctx.prec=precision;D=Decimal;n,ar,ai,co,si,ql,qr=constants()
        a=jacobian([ar]*N,[ai]*N,co,si);rhs=[[D(0)]*4 for _ in range(2*N)]
        for site,q,offset in [(0,ql,0),(N-1,qr,2)]:
            rhs[2*site][offset]=-q[0];rhs[2*site+1][offset]=-q[1]
            rhs[2*site][offset+1]=q[1];rhs[2*site+1][offset+1]=-q[0]
        solved=solve_real_banded(a,rhs)
        z=np.array([[float(v) for v in row] for row in solved]);response=z[::2]+1j*z[1::2]
        common=response[:,0]+response[:,2]
        return dict(N=N,precision=precision,density=float(n),alpha0=[float(ar),float(ai)],
                    response_real=response.real.tolist(),response_imag=response.imag.tolist(),
                    common_peak=float(max(abs(common))),common_relative_peak=float(max(abs(common)))/np.sqrt(float(n)))

def nonlinear_check(N,eta,precision=100):
    with localcontext() as ctx:
        ctx.prec=precision;D=Decimal;n,ar,ai,co,si,ql,qr=constants();eta=D(str(eta))
        x=[ar]*N;y=[ai]*N;U=D('-.0002');delta=D('.4998');gamma=D('.5');epsilon=D('38.51')
        def force():
            out=[]
            for j in range(N):
                h=delta+2*U*(x[j]*x[j]+y[j]*y[j]);fr=epsilon+h*x[j]+gamma*y[j];fi=h*y[j]-gamma*x[j]
                if j>0:fr+=co*x[j-1]-si*y[j-1];fi+=si*x[j-1]+co*y[j-1]
                if j<N-1:fr+=co*x[j+1]+si*y[j+1];fi+=-si*x[j+1]+co*y[j+1]
                if j in [0,N-1]:
                    q=ql if j==0 else qr;fr+=(1+eta)*q[0];fi+=(1+eta)*q[1]
                out.extend([fr,fi])
            return out
        residual=None
        for iteration in range(30):
            f=force();residual=max(abs(v) for v in f)
            if residual<D('1e-60'):break
            step=solve_real_banded(jacobian(x,y,co,si),[[-v] for v in f])
            for j in range(N):x[j]+=step[2*j][0];y[j]+=step[2*j+1][0]
        assert residual<D('1e-55'),residual
        difference=[((x[j]-ar)**2+(y[j]-ai)**2).sqrt() for j in range(N)]
        return dict(N=N,eta=float(eta),precision=precision,iterations=iteration,residual=str(residual),
                    relative_amplitude_change=float(max(difference)/n.sqrt()),
                    relative_density_change=float(max(abs((x[j]*x[j]+y[j]*y[j])/n-1) for j in range(N))))
