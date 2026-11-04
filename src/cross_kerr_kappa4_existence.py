"""kappa=4: prove and numerically check absence of intermediate bulk roots."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

DEST=Path('data/generated/compensation_n50_u_equal_v_delta0_kappa4')
DEST.mkdir(parents=True,exist_ok=True)
N=50;J=1.;U=V=-.0002;DELTA=0.;KAPPA=4.;GAMMA=KAPPA/2;BETA=2*(U+V)

def main():
    phi=np.linspace(0,np.pi,181)
    epsilon=np.unique(np.r_[np.linspace(0,200,401),np.geomspace(200,1e6,201)])
    density=np.full((len(phi),len(epsilon)),np.nan);counts=np.zeros(density.shape,dtype=np.int8)
    max_relative_residual=0.
    for i,p in enumerate(phi):
        c=DELTA+U+2*J*np.cos(p)
        for j,e in enumerate(epsilon):
            roots=np.roots([BETA*BETA,2*c*BETA,c*c+GAMMA*GAMMA,-e*e])
            real=[float(z.real) for z in roots if abs(z.imag)<1e-8*max(1,abs(z.real)) and z.real>=-1e-8]
            counts[i,j]=len(real)
            assert len(real)==1,(p,e,roots)
            n=max(0,real[0]);density[i,j]=n
            relative=abs(n*((c+BETA*n)**2+GAMMA*GAMMA)-e*e)/max(1,e*e)
            max_relative_residual=max(max_relative_residual,float(relative))
    cmax=DELTA+U+2*J;cmin=DELTA+U-2*J
    lower_derivative=GAMMA*GAMMA-max(abs(cmin),abs(cmax))**2/3
    critical_kappa=2*cmax/np.sqrt(3)
    assert lower_derivative>0
    summary=dict(parameters=dict(N=N,J=J,U=U,V=V,delta=DELTA,kappa=KAPPA),
                 phi_rows=len(phi),epsilon_values=len(epsilon),root_tests=int(counts.size),
                 epsilon_range=[float(epsilon.min()),float(epsilon.max())],max_positive_effective_detuning=cmax,
                 required_detuning=float(np.sqrt(3)*GAMMA),uniform_derivative_lower_bound=float(lower_derivative),
                 maximum_kappa_for_possible_three_roots=float(critical_kappa),
                 minimum_physical_root_count=int(counts.min()),maximum_physical_root_count=int(counts.max()),
                 intermediate_roots_found=0,max_relative_cubic_residual=max_relative_residual)
    np.savez_compressed(DEST/'bulk_root_check.npz',phi=phi,epsilon=epsilon,unique_density=density,physical_root_count=counts)
    (DEST/'bulk_root_check.json').write_text(json.dumps(summary,indent=2))
    plot(phi,epsilon,density,summary);report(summary)
    print(json.dumps(summary,indent=2),flush=True)

def plot(phi,epsilon,density,summary):
    fig,ax=plt.subplots(figsize=(10,6),layout='constrained');ax.set_facecolor('#e7e7e7')
    ax.set(xlim=(0,200),ylim=(0,np.pi),xlabel=r'Drive $\epsilon/J$',ylabel=r'Hopping phase $\phi$',
           title='Compensated intermediate-state existence | N=50, U=V=-0.0002, Delta=0, kappa=4')
    ax.set_yticks(np.arange(7)*np.pi/6,['0',r'$\pi/6$',r'$\pi/3$',r'$\pi/2$',r'$2\pi/3$',r'$5\pi/6$',r'$\pi$'])
    ax.text(.5,.5,'No separate intermediate homogeneous solution\nfor any drive or phase\n\n'+r'$\max_\phi c(\phi)=1.9998 < \sqrt{3}\kappa/2=3.4641$',
            transform=ax.transAxes,ha='center',va='center',fontsize=14,bbox=dict(facecolor='white',alpha=.9,edgecolor='none',pad=12))
    ax.legend(handles=[Patch(facecolor='#e7e7e7',label='No intermediate root; stability test not applicable')],loc='upper right',fontsize=9)
    ax.grid(alpha=.15);fig.savefig('figures/compensated_stability_map_u_equal_v_delta0_kappa4.png',dpi=200)
    fig.savefig('figures/compensated_stability_map_u_equal_v_delta0_kappa4.pdf')
    fig,ax=plt.subplots(figsize=(9,5),layout='constrained')
    for deg in [0,30,60,90,180]:
        i=int(deg);ax.loglog(epsilon[1:],density[i,1:],label=f'{deg} degrees')
    ax.set(xlabel=r'Drive $\epsilon/J$',ylabel='unique homogeneous-root density',title='Root-existence check, kappa=4 (stability not assessed)');ax.legend();ax.grid(alpha=.15,which='both')
    fig.savefig('figures/kappa4_unique_uniform_root.png',dpi=180)

def report(summary):
    text=f'''## Equal U=V and Delta=0 with kappa=4: no intermediate target at any drive

Keep N=50, J=1, U=V=-0.0002 and Delta=0, and increase kappa from 1 to **4**. The intermediate-root compensation-stability map is now empty: **there is no separate intermediate homogeneous bulk root at any epsilon>=0 or 0<=phi<=pi**. Enlarging the drive scan cannot restore that target. This is absence of the branch, not an unstable BdG eigenmode, a failed preparation, or a failed numerical stationary solver.

![kappa=4 map: no intermediate homogeneous target anywhere, shown over an enlarged drive interval](../figures/compensated_stability_map_u_equal_v_delta0_kappa4.png)

### Analytical check of the entire parameter space

With the same U n^2 and bond-counting conventions as above, define c(phi)=Delta+U+2J cos(phi)=-0.0002+2 cos(phi), beta=2(U+V)=-0.0008 and gamma=kappa/2=2. The cubic input-output relation is

`epsilon^2 = F(n) = n[(c+beta n)^2+gamma^2]`.

Three positive roots require c>sqrt(3)gamma. But max_phi c=1.9998, while sqrt(3)gamma={summary['required_detuning']:.10f}; the necessary condition fails at every phase. Equivalently,

`F'(n)=3(beta n+2c/3)^2+gamma^2-c^2/3`.

Since c ranges from -2.0002 to 1.9998, F'(n)>={summary['uniform_derivative_lower_bound']:.10f}>0 for all real n and all scanned phases. With F(0)=0 and F(n)->infinity as n->infinity, **exactly one nonnegative physical density root exists for every nonnegative drive**. At zero drive it is vacuum. Thus this conclusion is global in epsilon, not limited by the horizontal range shown in the plot.

The maximum loss rate for which a positive-density three-root region could exist with these fixed J, Delta and U is kappa_crit=2(Delta+U+2J)/sqrt(3)={summary['maximum_kappa_for_possible_three_roots']:.10f}. Kappa=4 exceeds it. V sets the nonlinear density/drive scales through beta, but increasing epsilon cannot overcome this detuning-versus-loss obstruction. Changing detuning or hopping would be a different parameter set and was not done.

### Numerical check over an enlarged drive range

As an independent numerical check, evaluate all roots of the cubic at 181 phases from 0 to pi, and {summary['epsilon_values']} drives: 0-200 on a linear grid of spacing 0.5, plus a logarithmic tail to **epsilon=1,000,000**. Across {summary['root_tests']} checks, the physical nonnegative-root count is always one; no intermediate target is found. The maximum relative cubic residual of the retained root is {summary['max_relative_cubic_residual']:.3g}. The plot displays 0-200 rather than the previous 0-85, and the analytic monotonicity argument covers drives beyond both the display and numerical check range.

![Unique homogeneous root versus drive over the extended numerical check](../figures/kappa4_unique_uniform_root.png)

These curves show the remaining single homogeneous root. Its root-specific endpoint compensation could be constructed using the cross-Kerr formula derived above, but **its BdG stability is not assessed in this intermediate-phase calculation**. No claim is made that all inhomogeneous mean-field solutions are absent, or that all compensated states are unstable. Only the distinct intermediate homogeneous branch is ruled out. Accordingly, no intermediate-state BdG spectra, adiabatic ramps, or covariance equations are calculated here.

### Reproduction and records

Run `python src/cross_kerr_kappa4_existence.py`. Parameters, analytic bounds and root-check statistics are in data/generated/compensation_n50_u_equal_v_delta0_kappa4/bulk_root_check.json. The NPZ companion stores the full phase/drive grids, physical root counts and unique-root densities. New figures use prefix compensated_stability_map_u_equal_v_delta0_kappa4; the kappa=1 maps above are retained for comparison.
'''
    p=Path('docs/uniform_topological_segment_stability.md');base=p.read_text(encoding='utf-8').split('## Equal U=V and Delta=0 with kappa=4: no intermediate target at any drive')[0].rstrip();p.write_text(base+'\n\n'+text,encoding='utf-8')

if __name__=='__main__':main()
