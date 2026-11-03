"""High-precision characteristic polynomials and Routh counts for banded real A.

Decimal calculations, not interval certificates. Caller controls precision.
"""
from decimal import Decimal, localcontext


def characteristic_banded(a, bandwidth=2):
    """Ascending coefficients of det(lambda I-A), by banded permutation DP."""
    n=len(a); w=bandwidth; zero=Decimal(0)
    states={0:[Decimal(1)]}
    for i in range(n):
        nxt={}
        for mask,p in states.items():
            for j in range(max(0,i-w),min(n,i+w+1)):
                bit=j-i+w
                if mask & (1<<bit): continue
                if j!=i and a[i][j]==0: continue
                updated=mask|(1<<bit)
                if i>=w and not updated&1: continue
                dest=updated>>1
                sign=Decimal(-1 if (mask>>(bit+1)).bit_count()%2 else 1)
                length=len(p)+(j==i)
                q=nxt.setdefault(dest,[zero]*length)
                if len(q)<length: q.extend([zero]*(length-len(q)))
                for k,v in enumerate(p):
                    q[k]-=sign*a[i][j]*v
                    if j==i: q[k+1]+=sign*v
        states=nxt
    return states[(1<<w)-1]


def count_right(coeff,shift):
    """Routh count of roots right of shift; positive row normalization."""
    zero=Decimal(0); c=Decimal(str(shift)); n=len(coeff)-1
    shifted=[coeff[-1]]
    for v in coeff[-2::-1]:
        new=[zero]*(len(shifted)+1)
        for j,z in enumerate(shifted): new[j]+=c*z; new[j+1]+=z
        new[0]+=v; shifted=new
    desc=shifted[::-1]
    a=desc[::2]; b=desc[1::2]
    b += [zero]*(len(a)-len(b))
    signs=[a[0]>0,b[0]>0]
    for _ in range(2,n+1):
        if not b[0]: raise ArithmeticError('Routh zero pivot: change shift/precision')
        row=[(b[0]*a[j+1]-a[0]*b[j+1])/b[0] for j in range(len(a)-1)]+[zero]
        scale=max(abs(v) for v in row)
        if not scale: raise ArithmeticError('Routh zero row: change shift/precision')
        row=[v/scale for v in row]
        signs.append(row[0]>0)
        a,b=b,row
    return sum(x!=y for x,y in zip(signs,signs[1:]))


def envelope_generator(n,width,kind):
    """Real, site-interleaved version of -i H in rotated Nambu coordinates."""
    one=Decimal(1); zero=Decimal(0)
    s=[]
    for j in range(n):
        d=Decimal(min(j+1,n-j))
        e=(2*d/Decimal(str(width))).exp() if width else None
        s.append((e-1)/(e+1) if width else one)
    a=[[zero]*(2*n) for _ in range(2*n)]
    for j in range(n):
        gs=s[j]**2 if kind=='paired' else one
        delta=Decimal('.1')*s[j]**2
        a[2*j][2*j]=gs; a[2*j+1][2*j+1]=-gs
        a[2*j][2*j+1]=delta; a[2*j+1][2*j]=-delta
        if j<n-1:
            gc=s[j]*s[j+1] if kind=='paired' else one
            a[2*j][2*j+2]=one+gc
            a[2*j+2][2*j]=gc-one
            a[2*j+1][2*j+3]=one-gc
            a[2*j+3][2*j+1]=-one-gc
    return a


def threshold(n,width,kind,precision=120,tolerance='1e-9'):
    with localcontext() as ctx:
        ctx.prec=precision
        a=envelope_generator(n,width,kind)
        coeff=characteristic_banded(a)
        lo=Decimal('1e-12'); hi=Decimal(4)
        if count_right(coeff,lo)==0:
            return dict(N=n,width=width,kind=kind,precision=precision,
                        growth_lower=0.,growth_upper=float(lo),critical_kappa=0.)
        assert count_right(coeff,hi)==0
        while hi-lo>Decimal(tolerance):
            mid=(lo+hi)/2
            if count_right(coeff,mid)>0: lo=mid
            else: hi=mid
        return dict(N=n,width=width,kind=kind,precision=precision,
                    growth_lower=float(lo),growth_upper=float(hi),critical_kappa=float(lo+hi))
