"""Fixed-background OBC comparison for locally topological mean-field segments."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eig
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from bulk_stability import matrix
from decimal import Decimal, localcontext
from banded_stability_decimal import characteristic_banded, count_right

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/generated/homogeneous_obc_n100'

def spectrum(n,d,g):
    h=matrix(n,1,np.pi/3,d,g,.5)
    w,vl,vr=eig(h,left=True,right=True)
    j=np.argmax(w.imag)
    cond=1/abs(np.vdot(vl[:,j],vr[:,j]))
    return dict(N=n,growth=float(w[j].imag),frequency=float(w[j].real),condition=float(cond),
                residual=float(np.linalg.norm(h@vr[:,j]-w[j]*vr[:,j])/np.linalg.norm(h)),
                sensitivity=float(np.finfo(float).eps*np.linalg.norm(h,np.inf)*cond))

def precise_count(n,d,g,precision):
    # A global Nambu rotation makes pairing imaginary without changing spectrum.
    # Real quadratures then have on-site [[g-gamma,d],[-d,-g-gamma]].
    with localcontext() as ctx:
        ctx.prec=precision
        D=Decimal; zero=D(0); half=D('.5'); s=D(3).sqrt()/2
        a=[[zero]*(2*n) for _ in range(2*n)]
        for j in range(n):
            a[2*j][2*j]=D(str(g))-half;a[2*j+1][2*j+1]=-D(str(g))-half
            a[2*j][2*j+1]=D(str(d));a[2*j+1][2*j]=-D(str(d))
            if j<n-1:
                for p,q,imag in [(j,j+1,-s),(j+1,j,s)]:
                    a[2*p][2*q]=imag;a[2*p+1][2*q+1]=imag
                    a[2*p][2*q+1]=half;a[2*p+1][2*q]=-half
        coeff=characteristic_banded(a,bandwidth=3)
        counts={str(x):count_right(coeff,D(str(x))) for x in [-.1,0,.1]}
        return dict(N=n,precision=precision,counts=counts)

def main():
    original=json.loads((OUT/'scan.json').read_text())[:161]
    refined=json.loads((OUT/'refinement.json').read_text())
    records={round(r['epsilon'],6):r for r in original if 35<r['epsilon']<45}
    records.update({round(r['epsilon'],6):r for r in refined})
    segments=[]
    for eps,r in sorted(records.items()):
        ids=np.where(np.array(r['nu'])!=0)[0]
        for group in np.split(ids,np.where(np.diff(ids)>1)[0]+1):
            if not len(group):continue
            n=np.array(r['density'])[group]
            row=dict(epsilon=eps,start=int(group[0]+1),end=int(group[-1]+1),length=len(group),
                     original_growth=r['growth'],density_min=float(n.min()),density_max=float(n.max()))
            row['tests']=[]
            for label,density in [('min',n.min()),('median',np.median(n)),('max',n.max())]:
                d=.4998-.0008*density;g=.0004*density
                for length in sorted(set([len(group),100])):
                    z=spectrum(length,d,g);z.update(sample=label,density=float(density),d=float(d),g=float(g))
                    # Real-k PBC maximum, continuous momentum (analytical here).
                    distance=max(abs(d)-1,0)
                    z['pbc_growth']=float(np.sqrt(max(g*g-distance*distance,0))-.5)
                    row['tests'].append(z)
            segments.append(row)
    (OUT/'uniform_segment_tests.json').write_text(json.dumps(segments,indent=2))
    mainseg=next(s for s in segments if s['epsilon']==38.5 and s['start']==4)
    checks=[]
    for sample in ['min','median','max']:
        z=next(z for z in mainseg['tests'] if z['sample']==sample and z['N']==100)
        for length in [50,100]:
            for precision in [100,160]:
                q=precise_count(length,z['d'],z['g'],precision);q['sample']=sample
                checks.append(q)
                print('precision check',q,flush=True)
                (OUT/'uniform_segment_precision.json').write_text(json.dumps(checks,indent=2))
    fig,ax=plt.subplots(2,1,figsize=(10,8),layout='constrained')
    for label in ['min','median','max']:
        rr=[(s,next(z for z in s['tests'] if z['sample']==label and z['N']==100)) for s in segments if s['length']>=5]
        ax[0].scatter([s['epsilon'] for s,z in rr],[z['growth'] for s,z in rr],s=12,label=label+' density')
    ax[0].axhline(0,color='k',ls='--');ax[0].legend();ax[0].set(xlabel='epsilon / J',ylabel='uniform N=100 OBC growth / J',title='Frozen coefficients from segments with at least 5 sites')
    for label in ['min','median','max']:
        z=next(z for z in mainseg['tests'] if z['sample']==label and z['N']==100)
        sizes=[5,10,20,30,50,75,100]
        ax[1].plot(sizes,[spectrum(n,z['d'],z['g'])['growth'] for n in sizes],'-o',label=label+' density')
    ax[1].axhline(0,color='k',ls='--');ax[1].legend();ax[1].set(xlabel='uniform OBC length',ylabel='growth / J',title='epsilon=38.50, segment 4-53')
    fig.savefig(ROOT/'figures/uniform_topological_segment_stability.png',dpi=180)
    print('main segment',mainseg,flush=True)

if __name__=='__main__':main()
