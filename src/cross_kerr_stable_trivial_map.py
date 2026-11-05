"""Full phi/drive map and resolved near-fold stable-trivial window, V=4U."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigvals
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator,FormatStrFormatter
from cross_kerr_trivial_search import operator,symbol,check
from verify_cross_kerr_trivial import precise

OUT=Path('data/generated/cross_kerr_trivial/map')
U=-.0002;V=4*U;DELTA=3.;GAMMA=.5;N=50

def background(phi,f):
    c=DELTA+U+2*np.cos(phi);rad=np.sqrt(c*c-3*GAMMA**2)
    t=(2*c-rad)/3+2*rad*f/3;n=t/(-2*(U+V))
    eps=np.sqrt(n*((c-t)**2+GAMMA**2))
    return n,eps

def growth(phi,f):
    n,_=background(phi,f)
    return float(eigvals(operator(N,n,phi,DELTA,1,V),check_finite=False).imag.max())

def grid(name,pe,fe):
    path=OUT/(name+'.npz')
    if path.exists():return dict(np.load(path))
    phi=(pe[:-1]+pe[1:])/2;f=(fe[:-1]+fe[1:])/2
    n,eps=background(phi[:,None],f[None,:]);rates=np.empty_like(n)
    a0,d0,b0=symbol(0,n,phi[:,None],DELTA,1,V)
    ap,dp,bp=symbol(np.pi,n,phi[:,None],DELTA,1,V)
    x0=(a0*a0-d0*d0+b0*b0).real;xp=(ap*ap-dp*dp+bp*bp).real
    assert np.all(x0>0)
    for i,p in enumerate(phi):
        rates[i]=[float(eigvals(operator(N,nn,p,DELTA,1,V),check_finite=False).imag.max()) for nn in n[i]]
        if i%20==0:print(name,'row',i,'/',len(phi),flush=True)
    _,xe=background(pe[:,None],fe[None,:]);ye=np.broadcast_to(pe[:,None],xe.shape)
    result=dict(phi=phi,fraction=f,density=n,epsilon=eps,growth=rates,X0=x0,Xpi=xp,xe=xe,ye=ye)
    np.savez_compressed(path,**result);return result

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    fe=np.unique(np.r_[0,np.geomspace(1e-7,.01,22),np.linspace(.02,.98,40),1-np.geomspace(1e-7,.01,22),1])
    full=grid('full',np.linspace(0,np.pi,121),fe)
    zoom=grid('zoom',np.deg2rad(np.linspace(34.9995,35.0005,81)),np.linspace(0,.004,65))
    stablechecks=[]
    for degrees in [34.9996,35,35.0004]:
        p=np.deg2rad(degrees);n,e=background(p,.001)
        r=dict(ratio=4,V=V,delta=DELTA,kappa=1,phi=float(p),density=float(n),epsilon=float(e),fraction=.001)
        r['check']=check(r,50);assert r['check']['growth'] < -1e-4
        r['precision']=[precise(r,prec) for prec in [100,160]]
        stablechecks.append(r)
    boundary=[]
    for p in np.deg2rad(np.linspace(34.9995,35.0005,41)):
        root=brentq(lambda f:growth(p,f),.001,.004,xtol=1e-12)
        _,eps=background(p,root);_,fold=background(p,0)
        boundary.append(dict(phi=float(p),fraction=float(root),epsilon=float(eps),fold=float(fold),drive_width=float(fold-eps)))
    summary=dict(parameters=dict(N=N,U=U,V=V,delta=DELTA,kappa=1,J=1,omega_ref=0),
                 checks=stablechecks,boundary=boundary,counts={})
    colors=['#26866d','#d77c54','#607da6','#a542a5','#edc64b']
    fig,axes=plt.subplots(1,2,figsize=(14,6));fig.subplots_adjust(left=.07,right=.96,bottom=.24,top=.80,wspace=.28)
    for ax,data,name in zip(axes,[full,zoom],['full','zoom']):
        rates=data['growth'];trivial=data['Xpi']>0;stable=rates < -1e-7
        label=np.where(stable,np.where(trivial,3,0),np.where(trivial,2,1));label[abs(rates)<=1e-7]=4
        summary['counts'][name]={str(i):int((label==i).sum()) for i in range(5)}
        ax.set_facecolor('#ededed')
        ax.pcolormesh(data['xe'],np.rad2deg(data['ye']),label,cmap=ListedColormap(colors),norm=BoundaryNorm(np.arange(6)-.5,5),shading='flat',rasterized=True)
        ax.plot(data['xe'][:,0],np.rad2deg(data['ye'][:,0]),'k',lw=1)
        if name=='full':
            ax.plot(data['xe'][:,-1],np.rad2deg(data['ye'][:,-1]),'k',lw=1)
            ax.contour(data['epsilon'],np.broadcast_to(np.rad2deg(data['phi'])[:,None],rates.shape),data['Xpi'],levels=[0],colors='#333',linestyles='--',linewidths=.7)
            ax.plot(87.095716964,35,'*',color='#f7ce47',mec='k',ms=12,zorder=5)
            ax.annotate('Zoom near this point',xy=(87.095716964,35),xytext=(92,65),arrowprops=dict(arrowstyle='->'),fontsize=10)
            ax.set(xlim=(0,125),ylim=(0,180),title='Full intermediate-branch diagram')
        else:
            # The computed near-fold band extends below the displayed drive
            # window; the upper boundary is the exact root-existence fold.
            by=np.array([np.rad2deg(b['phi']) for b in boundary]);bx=[b['epsilon'] for b in boundary]
            folds=[b['fold'] for b in boundary]
            ax.fill_betweenx(by,87.0948,bx,color=colors[2],zorder=2)
            ax.fill_betweenx(by,bx,folds,color=colors[3],zorder=2)
            ax.fill_betweenx(by,folds,87.0965,color='#ededed',zorder=2)
            ax.plot(folds,by,color='k',lw=1,zorder=3)
            ax.plot(bx,by,color='#222',ls=':',lw=1.3,label='Zero-growth boundary')
            ax.plot(87.095716964,35,'*',color='#f7ce47',mec='k',ms=12,zorder=5)
            ax.set(xlim=(87.0948,87.0965),ylim=(34.9995,35.0005),title=r'Resolved stable trivial window ($\nu=0$)')
            ax.xaxis.set_major_locator(MaxNLocator(4));ax.xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
            ax.yaxis.set_major_locator(MaxNLocator(5));ax.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))
        ax.set(xlabel=r'Drive $\epsilon/J$',ylabel=r'Hopping phase $\phi$ (degrees)');ax.grid(alpha=.12)
    handles=[Patch(color=colors[0],label=r'Stable, $\nu=-1$'),Patch(color=colors[1],label=r'Unstable, $\nu=-1$'),Patch(color=colors[2],label=r'Unstable, $\nu=0$'),Patch(color=colors[3],label=r'Stable, $\nu=0$'),Patch(color='#ededed',label='No intermediate root')]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.10),ncol=5,fontsize=10,frameon=False)
    fig.suptitle(r'Stable trivial compensated intermediate states with cross-Kerr coupling',fontsize=15,y=.97)
    fig.text(.5,.89,r'$N=50$, $J=\kappa=1$, $\Delta=3$, $U=-2\times10^{-4}$, $V=4U$; zero-frequency bulk winding',ha='center',fontsize=12)
    fig.text(.5,.04,'Purple: stable and topologically trivial. Black solid: branch merger (gapless). Dotted: OBC stability boundary.',ha='center',fontsize=10)
    for ext in ['png','svg']:fig.savefig(f'figures/cross_kerr_stable_trivial_phase_map.{ext}',dpi=200)
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(dict(counts=summary['counts'],boundary_example=boundary[len(boundary)//2]),indent=2))

if __name__=='__main__':main()
