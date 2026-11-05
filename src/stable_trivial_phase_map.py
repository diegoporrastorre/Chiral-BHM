"""Joint OBC stability/topology map and no-stable-trivial verification, V=0."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json
import numpy as np
from scipy.linalg import eigvals
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib.patches import Patch
from bulk_stability import matrix

def main():
    dest=Path('data/generated/stable_trivial');dest.mkdir(parents=True,exist_ok=True)
    N=50;a=3.;gamma=.5;U=-.0002
    pe=np.linspace(0,np.pi,161);phi=(pe[:-1]+pe[1:])/2
    fe=np.sin(np.linspace(0,np.pi/2,65))**2;f=(fe[:-1]+fe[1:])/2
    c=a+2*np.cos(phi);rad=np.sqrt(c*c-3*gamma**2)
    t=(2*c[:,None]-rad[:,None])/3+2*rad[:,None]*f/3
    n=t/(-2*U);eps=np.sqrt(n*((c[:,None]-t)**2+gamma**2));d=a-2*t
    x0=t*t-gamma**2-(d+2*np.cos(phi[:,None]))**2
    xp=t*t-gamma**2-(d-2*np.cos(phi[:,None]))**2
    growth=np.zeros_like(n);counts=np.zeros_like(n,dtype=int)
    path=dest/'phase_map.npz'
    if path.exists():
        saved=np.load(path);growth=saved['growth'];counts=saved['growing_modes']
    else:
        for i,p in enumerate(phi):
            for j in range(len(f)):
                w=eigvals(matrix(N,1,p,d[i,j],t[i,j],gamma),check_finite=False)
                growth[i,j]=w.imag.max();counts[i,j]=(w.imag>1e-7).sum()
            if i%20==0:print('map row',i,'/160',flush=True)
        np.savez_compressed(path,phi=phi,density=n,epsilon=eps,growth=growth,growing_modes=counts,X0=x0,Xpi=xp)
    trivial=xp>0;stable=growth < -1e-7
    assert np.all(x0>0)
    assert np.all(counts[trivial]==N)
    # Independent matrix-basis check at representative trivial backgrounds,
    # and selected sizes, with exact norm inequality verified.
    checks=[]
    for i,j in [(20,20),(80,32),(130,20)]:
        if not trivial[i,j]:continue
        for size in [8,20,50,80]:
            p=phi[i];dd=d[i,j];g=t[i,j]
            D=dd*np.eye(size)+np.cos(p)*(np.eye(size,k=1)+np.eye(size,k=-1))
            K=np.sin(p)*(np.eye(size,k=-1)-np.eye(size,k=1))
            G=np.block([[K+(g-gamma)*np.eye(size),D],[-D,K-(g+gamma)*np.eye(size)]])
            vals=eigvals(G.astype(complex));norm=float(np.linalg.norm(D,2))
            assert norm**2<g*g-gamma*gamma
            assert (vals.real>1e-7).sum()==size
            ref=eigvals(matrix(size,1,p,dd,g,gamma))
            assert abs(vals.real.max()-ref.imag.max())<1e-8
            checks.append(dict(N=size,phi=float(p),density=float(n[i,j]),growing_modes=size,norm_margin=float(g*g-gamma*gamma-norm*norm),growth=float(vals.real.max())))
    summary=dict(parameters=dict(N=N,J=1,U=U,V=0,delta=3.0002,kappa=1),samples=n.size,
                 stable_topological=int((stable & ~trivial).sum()),unstable_topological=int(((growth>1e-7)&~trivial).sum()),
                 stable_trivial=int((stable&trivial).sum()),unstable_trivial=int(((growth>1e-7)&trivial).sum()),checks=checks)
    (dest/'phase_map.json').write_text(json.dumps(summary,indent=2))
    label=np.where(stable,np.where(trivial,3,0),np.where(trivial,2,1))
    label[abs(growth)<=1e-7]=4
    ce=a+2*np.cos(pe);re=np.sqrt(ce*ce-3*gamma**2)
    te=(2*ce[:,None]-re[:,None])/3+2*re[:,None]*fe/3
    xe=np.sqrt(te/(-2*U)*((ce[:,None]-te)**2+gamma**2));ye=np.broadcast_to(pe[:,None],xe.shape)
    colors=['#26866d','#cf6b4b','#547ca8','#b051b0','#edc64b']
    fig,ax=plt.subplots(figsize=(10,7));fig.subplots_adjust(left=.10,right=.97,bottom=.12,top=.84)
    ax.set_facecolor('#e7e7e7')
    ax.pcolormesh(xe,ye,label,cmap=ListedColormap(colors),norm=BoundaryNorm(np.arange(6)-.5,5),shading='flat',rasterized=True)
    ax.plot(xe[:,0],pe,'k',lw=.8);ax.plot(xe[:,-1],pe,'k',lw=.8)
    ax.contour(eps,np.broadcast_to(phi[:,None],eps.shape),xp,levels=[0],colors='#472257',linewidths=1,linestyles='--')
    ax.set(xlabel=r'Drive $\epsilon/J$',ylabel=r'Hopping phase $\phi$',xlim=(0,220),ylim=(0,np.pi))
    ax.set_yticks(np.arange(7)*np.pi/6,['0',r'$\pi/6$',r'$\pi/3$',r'$\pi/2$',r'$2\pi/3$',r'$5\pi/6$',r'$\pi$'])
    ax.grid(alpha=.12)
    ax.legend(handles=[Patch(color=colors[0],label=r'Stable, topological ($\nu=-1$)'),Patch(color=colors[1],label=r'Unstable, topological ($\nu=-1$)'),Patch(color=colors[2],label=r'Unstable, trivial ($\nu=0$)'),Patch(color='#e7e7e7',label='No distinct intermediate root')],loc='upper right',fontsize=9)
    fig.suptitle('Compensated intermediate branch: no stable trivial region',fontsize=15,y=.97)
    ax.set_title(r'$N=50$, $J=\kappa=1$, $U=-2\times10^{-4}$, $V=0$, $\Delta=3.0002$; $\omega_0=0$',fontsize=11,pad=15)
    ax.text(.035,.06,'Stable + trivial: excluded analytically for this model.\nDashed curves: point-gap closure at '+r'$k=\pi$'+'.',transform=ax.transAxes,fontsize=10,bbox=dict(facecolor='white',edgecolor='none',alpha=.95))
    for ext in ['png','svg']:fig.savefig(f'figures/compensated_stability_topology_detuning3.{ext}',dpi=200)
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
