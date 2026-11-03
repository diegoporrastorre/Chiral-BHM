"""N=100: sharp left edge, exactly uniform core, right-only boundary layer."""
from pathlib import Path
from decimal import Decimal,localcontext
import json,time
from banded_stability_decimal import characteristic_banded,count_right

ROOT=Path(__file__).resolve().parents[1]
PATH=ROOT/'data/generated/bdg_right_absorber.json'


def profiles(n,length,family,eta=0):
    one=Decimal(1); zero=Decimal(0)
    s=[]; loss=[]
    for j in range(n):
        x=max(zero,Decimal(j-(n-length)+1)/Decimal(length)) if length else zero
        if family.startswith('cap'):
            power=int(family[3:]); f=1-x**power
        elif family in ['quintic','loss']:
            f=1-10*x**3+15*x**4-6*x**5
        elif family=='cubic': f=1-3*x**2+2*x**3
        elif family.startswith('tanh'):
            shape=Decimal(family[4:])
            e=(2*shape*(1-x)).exp(); e0=(2*shape).exp()
            f=(e-1)/(e+1)*(e0+1)/(e0-1)
        elif family=='linear': f=1-x
        else: raise ValueError(family)
        s.append(one if family=='loss' or family.startswith('cap') else f)
        loss.append(Decimal(str(eta))*(1-f))
    return s,loss


def generator(n,length,family,eta=0):
    s,loss=profiles(n,length,family,eta)
    zero=Decimal(0)
    a=[[zero]*(2*n) for _ in range(2*n)]
    for j in range(n):
        g=s[j]**2; delta=Decimal('.1')*g
        a[2*j][2*j]=g-loss[j]; a[2*j+1][2*j+1]=-g-loss[j]
        a[2*j][2*j+1]=delta; a[2*j+1][2*j]=-delta
        if j<n-1:
            gc=s[j]*s[j+1]
            a[2*j][2*j+2]=1+gc; a[2*j+2][2*j]=gc-1
            a[2*j+1][2*j+3]=1-gc; a[2*j+3][2*j+1]=-1-gc
    return a


def spectral_abscissa(n,length,family,eta=0,precision=160,tol='1e-5'):
    with localcontext() as ctx:
        ctx.prec=precision
        coeff=characteristic_banded(generator(n,length,family,eta))
        def val(x):
            y=Decimal(0)
            for c in coeff[::-1]: y=y*x+c
            return y
        hi=Decimal(4); lo=hi
        while lo>Decimal('.01') and val(lo)>=0:
            hi=lo; lo-=Decimal('.01')
        if val(lo)<0:
            while hi-lo>Decimal(tol):
                mid=(lo+hi)/2
                if val(mid)>0: hi=mid
                else: lo=mid
            counts=[count_right(coeff,lo),count_right(coeff,hi)]
        else: counts=[0,1]
        if counts[0]==0 or counts[1]!=0:
            lo=Decimal('1e-9'); hi=Decimal(4)
            assert count_right(coeff,hi)==0
            if count_right(coeff,lo)==0:
                return dict(N=n,length=length,family=family,eta=eta,critical_kappa=0.,
                            lower=0.,upper=float(lo),precision=precision,counts=[0,0])
            while hi-lo>Decimal(tol):
                mid=(lo+hi)/2
                if count_right(coeff,mid)>0: lo=mid
                else: hi=mid
            counts=[count_right(coeff,lo),count_right(coeff,hi)]
        return dict(N=n,length=length,family=family,eta=eta,critical_kappa=float(lo+hi),
                    lower=float(lo),upper=float(hi),precision=precision,counts=counts)


def verify(row,precision=260):
    with localcontext() as ctx:
        ctx.prec=precision
        coeff=characteristic_banded(generator(row['N'],row['length'],row['family'],row['eta']))
        counts=[count_right(coeff,Decimal(str(row[k]))+Decimal(sign)*Decimal('1e-10'))
                for k,sign in [('lower',-1),('upper',1)]]
        assert counts[0]>0 and counts[1]==0,(row,counts)
    row.update(check_precision=precision,check_counts=counts)
    return row


def run():
    result=json.loads(PATH.read_text()) if PATH.exists() else dict(parameters=dict(N=100,gs=1,gc=1,Delta=.1,J=1,phi='pi/2',P=0),rows=[])
    cases=[(0,'quintic',0)]+[(l,f,0) for l in [20,40,60,80,90] for f in ['linear','cubic','quintic','tanh2','tanh4']]
    cases += [(l,'loss',e) for l in [40,80] for e in [1,3,10]]
    cases += [(l,'quintic',0) for l in [92,94,96,98,99]]
    cases += [(l,f,e) for l in [80,90] for f in ['cap1','cap2','cap4','cap8'] for e in [1,3,10]]
    cases += [(l,'quintic',e) for l in [60,80,90] for e in [.3,1,3,10]]
    for length,family,eta in cases:
        if any(r['length']==length and r['family']==family and r['eta']==eta for r in result['rows']): continue
        start=time.perf_counter()
        row=verify(spectral_abscissa(100,length,family,eta))
        row['elapsed_seconds']=time.perf_counter()-start
        result['rows'].append(row)
        PATH.write_text(json.dumps(result,indent=2))
        print('DONE',length,family,eta,row['critical_kappa'],row['elapsed_seconds'],flush=True)
    return result


if __name__=='__main__': run()
