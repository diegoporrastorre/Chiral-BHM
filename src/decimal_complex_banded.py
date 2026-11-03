"""Small Decimal-complex implementation for precision-refined banded solves."""
from decimal import Decimal


class C:
    def __init__(self,re=0,im=0):
        if isinstance(re,C): self.re,self.im=re.re,re.im
        else: self.re,self.im=Decimal(str(re)),Decimal(str(im))
    def __add__(self,b):
        b=C(b); return C(self.re+b.re,self.im+b.im)
    __radd__=__add__
    def __neg__(self): return C(-self.re,-self.im)
    def __sub__(self,b): return self+-C(b)
    def __rsub__(self,b): return C(b)+-self
    def __mul__(self,b):
        b=C(b); return C(self.re*b.re-self.im*b.im,self.re*b.im+self.im*b.re)
    __rmul__=__mul__
    def norm2(self): return self.re*self.re+self.im*self.im
    def __truediv__(self,b):
        b=C(b); d=b.norm2()
        return C((self.re*b.re+self.im*b.im)/d,(self.im*b.re-self.re*b.im)/d)
    def __rtruediv__(self,b): return C(b)/self
    def sqrt(self):
        r=self.norm2().sqrt()
        x=max(Decimal(0),(r+self.re)/2).sqrt()
        y=max(Decimal(0),(r-self.re)/2).sqrt()
        return C(x,-y if self.im<0 else y)
    def magnitude(self): return self.norm2().sqrt()


def solve_banded(a,b,width=2):
    """Partial-pivot Gaussian elimination: lower width w, upper fill <=2w."""
    n=len(b)
    a=[list(row) for row in a]; b=list(b)
    for k in range(n):
        pivot=max(range(k,min(n,k+width+1)),key=lambda i:a[i][k].norm2())
        if a[pivot][k].norm2()==0: raise ArithmeticError('Singular banded system')
        if pivot!=k: a[k],a[pivot]=a[pivot],a[k]; b[k],b[pivot]=b[pivot],b[k]
        for i in range(k+1,min(n,k+width+1)):
            factor=a[i][k]/a[k][k]
            for j in range(k+1,min(n,k+2*width+1)):
                a[i][j]=a[i][j]-factor*a[k][j]
            b[i]=b[i]-factor*b[k]
            a[i][k]=C(0)
    x=[C(0) for _ in range(n)]
    for i in range(n-1,-1,-1):
        rhs=b[i]
        for j in range(i+1,min(n,i+2*width+1)): rhs=rhs-a[i][j]*x[j]
        x[i]=rhs/a[i][i]
    return x
