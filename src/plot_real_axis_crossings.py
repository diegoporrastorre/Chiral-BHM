"""Illustrate determinant crossings on one intermediate branch, V=0."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def main():
    U=-.0002; delta=3.0002; phi=np.pi/3; gamma=.5
    a=delta+U;h=2*np.cos(phi);c=a+h
    transition=(8+np.sqrt(13))/6
    ts=[2.3,transition,1.6]
    titles=[r'Topological: $\nu(0)=-1$',r'Transition: $\nu(0)$ undefined',r'Trivial: $\nu(0)=0$']
    fig,axes=plt.subplots(1,3,figsize=(15,5.7),sharex=True,sharey=True)
    fig.subplots_adjust(left=.055,right=.985,bottom=.29,top=.73,wspace=.13)
    results=[]
    for index,(ax,t,title) in enumerate(zip(axes,ts,titles)):
        n=t/(-2*U);d=a-2*t;eps=np.sqrt(n*((c-t)**2+gamma**2))
        roots=np.roots([4*U**2,4*c*U,c*c+gamma*gamma,-eps**2])
        assert np.max(abs(roots.imag))<1e-8
        assert abs(np.sort(roots.real)[1]-n)<1e-6
        def det(k):
            return (2*np.sin(phi)*np.sin(k)-1j*gamma)**2-(d+h*np.cos(k))**2+t*t
        x0=float(det(0).real);xp=float(det(np.pi).real)
        assert x0>0
        if index==1: assert abs(xp)<1e-12
        else:
            k=np.linspace(-np.pi,np.pi,16384,endpoint=False);v=det(k)
            nu=int(round(np.angle(np.roll(v,-1)/v).sum()/(2*np.pi)))
            assert nu==(-1 if index==0 else 0)
        for lo,hi,col in [(-np.pi,0,'#2475ad'),(0,np.pi,'#d36a2d')]:
            k=np.linspace(lo,hi,2001);z=det(k)
            ax.plot(z.real,z.imag,color=col,lw=2)
            # Arrows follow increasing momentum, not increasing real coordinate.
            for midpoint in [lo+(hi-lo)*.32,lo+(hi-lo)*.72]:
                z1,z2=det(midpoint-.035),det(midpoint+.035)
                ax.annotate('',xy=(z2.real,z2.imag),xytext=(z1.real,z1.imag),
                            arrowprops=dict(arrowstyle='->',color=col,lw=1.8,mutation_scale=13))
        ax.axhline(0,color='#888',lw=.7);ax.axvline(0,color='#888',lw=.7)
        ax.plot(0,0,'x',color='#c32469',ms=9,mew=2,zorder=5)
        ax.plot(x0,0,'o',color='#222',ms=6,zorder=6)
        ax.plot(xp,0,'s',color='#7244a0',ms=6,zorder=6)
        ax.annotate(r'$X_0$'+f' = {x0:.3f}',(x0,0),xytext=(0,23),textcoords='offset points',ha='center',fontsize=10,
                    arrowprops=dict(arrowstyle='-',color='#555'))
        ax.annotate(r'$X_\pi$'+f' = {xp:.3f}',(xp,0),xytext=(0,-35),textcoords='offset points',ha='center',fontsize=10,
                    arrowprops=dict(arrowstyle='-',color='#7244a0'))
        ax.set_title(title+'\n'+rf'$n={n:.2f},\quad\epsilon/J={eps:.4f}$',fontsize=11,pad=15)
        ax.set(xlim=(-2.7,7),ylim=(-2.35,2.35),xlabel=r'$\mathrm{Re}\,D(k)/J^2$')
        ax.set_aspect('equal',adjustable='box');ax.grid(alpha=.12)
        results.append(dict(density=n,epsilon=float(eps),t=float(t),X0=x0,Xpi=xp,winding=None if index==1 else nu))
    axes[0].set_ylabel(r'$\mathrm{Im}\,D(k)/J^2$')
    fig.suptitle(r'Real-axis crossings of $D(k)=\det[-\mathcal{H}(k)]$ at $\omega_0=0$',fontsize=16,y=.97)
    fig.text(.5,.89,r'One homogeneous intermediate branch: $J=\kappa=1$, $U=-2\times10^{-4}$, $V=0$, $\Delta=3.0002$, $\phi=\pi/3$',ha='center',fontsize=11)
    fig.legend(handles=[Line2D([0],[0],color='#2475ad',lw=2,label=r'$-\pi\to0$'),Line2D([0],[0],color='#d36a2d',lw=2,label=r'$0\to\pi$'),Line2D([0],[0],marker='o',color='none',markerfacecolor='#222',label=r'$k=0$'),Line2D([0],[0],marker='s',color='none',markerfacecolor='#7244a0',label=r'$k=\pm\pi$'),Line2D([0],[0],marker='x',color='#c32469',ls='none',label='Origin')],loc='lower center',bbox_to_anchor=(.5,.13),ncol=5,frameon=False)
    fig.text(.5,.075,r'$X_0>0$ on every intermediate state.  $X_\pi<0$: encloses the origin;  $X_\pi=0$: gap closes;  $X_\pi>0$: no winding.',ha='center',fontsize=11)
    for extension in ['png','svg']:
        fig.savefig(f'figures/intermediate_real_axis_crossings.{extension}',dpi=200)
    dest=Path('data/generated/intermediate_topology_general');dest.mkdir(parents=True,exist_ok=True)
    (dest/'real_axis_crossings.json').write_text(json.dumps(results,indent=2))
    print(json.dumps(results,indent=2))


if __name__=='__main__':main()

