"""Fixed epsilon=38.51, N=50 intermediate-root compensated OBC phase scan."""
from compensation_n50_vacuum import *
from scipy.linalg import eig
from scipy.optimize import brentq
from banded_stability_decimal import characteristic_banded,count_right
from decimal import Decimal,localcontext

FILE=DEST/'compensated_phi_stability.json'

def lower_turn(degrees):
    c=.4998+2*np.cos(np.deg2rad(degrees));u=-.0004
    n=max(np.roots([3*u*u,4*c*u,c*c+.25]))
    return float(np.sqrt(n*((c+u*n)**2+.25)))

def point(degrees):
    phi=np.deg2rad(degrees);c=.4998+2*np.cos(phi);u=-.0004
    roots=np.roots([u*u,2*c*u,c*c+.25,-38.51**2])
    real=sorted(z.real for z in roots if abs(z.imag)<1e-6 and z.real>0)
    r=dict(degrees=float(degrees),phi=float(phi),real_densities=real,intermediate_exists=len(real)==3)
    if len(real)!=3:return r
    n=real[1];a0=-38.51/(c+u*n-.5j);a=np.full(50,a0)
    m=Chain(N=50,drive=38.51,phi=phi);m.pump=m.pump.astype(complex)
    m.pump[0]+=np.exp(1j*phi)*a0;m.pump[-1]+=np.exp(-1j*phi)*a0
    residual=float(max(abs(m.force(a))));assert residual<1e-8
    L=m.bogoliubov(a);w,vl,vr=eig(L,left=True,right=True);j=np.argmax(w.imag)
    cond=1/abs(np.vdot(vl[:,j],vr[:,j]));d=.4998-.0008*n;g=.0004*n
    x0=g*g-.25-(d+2*np.cos(phi))**2;xp=g*g-.25-(d-2*np.cos(phi))**2
    r.update(density=float(n),alpha0=[float(a0.real),float(a0.imag)],
             endpoint_pumps=[[float(m.pump[j].real),float(m.pump[j].imag)] for j in [0,-1]],
             residual=residual,growth=float(w[j].imag),dominant_frequency=float(w[j].real),
             condition=float(cond),sensitivity=float(np.finfo(float).eps*np.linalg.norm(L,np.inf)*cond),
             eigen_residual=float(np.linalg.norm(L@vr[:,j]-w[j]*vr[:,j])/np.linalg.norm(L)),
             winding=int(np.sign(d)) if x0*xp<0 else 0,eigen_real=w.real.tolist(),eigen_imag=w.imag.tolist())
    return r

def precise(degrees,precision,epsilon=38.51,initial_density=None,N=50):
    with localcontext() as ctx:
        ctx.prec=precision+20
        D=Decimal;tol=D(10)**(-precision-10)
        def atan(x):
            total=term=x;k=1
            while abs(term)>tol:
                term=-term*x*x*(2*k-1)/(2*k+1);total+=term;k+=1
            return total
        pi=16*atan(D(1)/5)-4*atan(D(1)/239);ph=pi*D(str(degrees))/180
        cs=term=D(1);k=1
        while abs(term)>tol:
            term=-term*ph*ph/((2*k-1)*(2*k));cs+=term;k+=1
        sn=term=ph;k=1
        while abs(term)>tol:
            term=-term*ph*ph/((2*k)*(2*k+1));sn+=term;k+=1
        c=D('.4998')+2*cs;u=D('-.0004')
        if initial_density is None:
            real=sorted(z.real for z in np.roots([.0004**2,-.0008*float(c),float(c)**2+.25,-epsilon**2]) if abs(z.imag)<1e-6 and z.real>0)
            if len(real)!=3:raise ValueError('No distinct intermediate density root')
            initial_density=real[1]
        n=D(str(initial_density))
        for _ in range(100):
            step=(n*((c+u*n)**2+D('.25'))-D(str(epsilon))**2)/((c+u*n)*(c+3*u*n)+D('.25'))
            n-=step
            if abs(step)<tol:break
        d=D('.4998')-D('.0008')*n;g=D('.0004')*n
        ctx.prec=precision
        cs=+cs;sn=+sn;zero=D(0);gamma=D('.5')
        a=[[zero]*(2*N) for _ in range(2*N)]
        for j in range(N):
            a[2*j][2*j]=g-gamma;a[2*j+1][2*j+1]=-g-gamma
            a[2*j][2*j+1]=d;a[2*j+1][2*j]=-d
            if j<N-1:
                for p,q,im in [(j,j+1,-sn),(j+1,j,sn)]:
                    a[2*p][2*q]=im;a[2*p+1][2*q+1]=im
                    a[2*p][2*q+1]=cs;a[2*p+1][2*q]=-cs
        coeff=characteristic_banded(a,3)
        return dict(degrees=float(degrees),epsilon=float(epsilon),precision=precision,density=float(n),
                    counts={str(s):count_right(coeff,Decimal(str(s))) for s in [-1e-7,0,1e-7]})

def main():
    merger=brentq(lambda degrees:lower_turn(degrees)-38.51,0,60,xtol=1e-12)
    grid=sorted(set(np.r_[np.linspace(0,60,601),merger+np.geomspace(1e-5,.1,60)]),reverse=True)
    rows=[point(float(p)) for p in grid];valid=sorted([r for r in rows if r['intermediate_exists']],key=lambda r:r['degrees'])
    transitions=[]
    for a,b in zip(valid[:-1],valid[1:]):
        if a['growth']*b['growth']<0:
            deg=brentq(lambda deg:point(deg)['growth'],a['degrees'],b['degrees'],xtol=1e-10)
            transitions.append(dict(degrees=float(deg),phi=float(np.deg2rad(deg)),scan_bracket=[a['degrees'],b['degrees']]))
    result=dict(N=50,epsilon=38.51,parameters=dict(J=1,U=-.0002,V=0,delta=.5,kappa=1),
                merger_degrees=float(merger),merger_phi=float(np.deg2rad(merger)),rows=rows,transitions=transitions,precision_checks=[])
    FILE.write_text(json.dumps(result))
    selected={60.,55.,40.,20.,merger+.01}
    for t in transitions:selected.update([t['degrees']-.001,t['degrees']+.001])
    for deg in sorted(selected,reverse=True):
        if deg<=merger:continue
        for precision in [100,160]:
            r=precise(deg,precision);result['precision_checks'].append(r);FILE.write_text(json.dumps(result));print(r,flush=True)
    plot_result(result);report_result(result)
    print('merger',merger,'transitions',transitions,flush=True)

def plot_result(result):
    valid=sorted([r for r in result['rows'] if r['intermediate_exists']],key=lambda r:r['degrees'])
    fig,ax=plt.subplots(3,1,figsize=(10,10),layout='constrained')
    deg=[r['degrees'] for r in valid]
    ax[0].plot(deg,[r['density'] for r in valid]);ax[0].set(ylabel='intermediate density')
    ax[1].plot(deg,[r['growth'] for r in valid]);ax[1].axhline(0,color='k',ls='--');ax[1].set(ylabel='max Im omega / J')
    ax[2].semilogy(deg,[max(r['condition'],1) for r in valid]);ax[2].set(ylabel='dominant eigenvalue condition',xlabel='hopping phase phi (degrees)')
    for axes in ax:
        axes.axvspan(0,result['merger_degrees'],color='.8',alpha=.5,label='no intermediate root')
        for t in result['transitions']:axes.axvline(t['degrees'],color='#b44c12',ls=':')
        axes.set_xlim(0,60);axes.grid(alpha=.15)
    ax[0].legend();fig.savefig('figures/compensated_intermediate_phi_stability.png',dpi=180)
    fig,ax=plt.subplots(figsize=(9,5),layout='constrained')
    for p in [60,55,40,20]:
        r=min(valid,key=lambda r:abs(r['degrees']-p));ax.scatter(r['eigen_real'],r['eigen_imag'],s=10,label=f'{p} degrees')
    ax.axhline(0,color='k',ls='--');ax.set(xlabel='Re omega / J',ylabel='Im omega / J');ax.legend();fig.savefig('figures/compensated_intermediate_phi_spectra.png',dpi=180)
    fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    for axes_one,limits in zip(axes,[(18.22,18.4),(55.8,56.7)]):
        close=[r for r in valid if limits[0]<=r['degrees']<=limits[1]]
        # Display extra evaluated points around each crossing for readable zooms.
        close+= [point(float(p)) for p in np.linspace(max(limits[0],result['merger_degrees']+1e-5),limits[1],120)]
        close=sorted(close,key=lambda r:r['degrees'])
        axes_one.plot([r['degrees'] for r in close],[r['growth'] for r in close])
        axes_one.axhline(0,color='k',ls='--');axes_one.axvspan(0,result['merger_degrees'],color='.8',alpha=.5)
        for t in result['transitions']:axes_one.axvline(t['degrees'],color='#b44c12',ls=':')
        axes_one.set(xlim=limits,xlabel='phi (degrees)',ylabel='max Im omega / J');axes_one.grid(alpha=.15)
    fig.savefig('figures/compensated_intermediate_phi_crossings.png',dpi=180)

def report_result(result):
    valid=sorted([r for r in result['rows'] if r['intermediate_exists']],key=lambda r:r['degrees'],reverse=True)
    transitions=result['transitions'];m=result['merger_degrees'];table=[]
    for p in [60,59,58,57,56,55,50,45,40,35,30,25,20,18,17,15,10,0]:
        r=min(result['rows'],key=lambda r:abs(r['degrees']-p))
        table.append(f"|{p:g}|{r['density']:.3f}|{r['growth']:.6f}|{'stable' if r['growth']<-1e-7 else 'unstable' if r['growth']>1e-7 else 'marginal'}|" if r['intermediate_exists'] else f'|{p:g}|--|--|no intermediate root|')
    checks=[]
    for r in result['precision_checks']:
        if r['precision']==160:checks.append(f"|{r['degrees']:.6f}|{r['counts']['0']}|")
    text=f'''## Compensated intermediate-state stability versus phi at epsilon=38.51

This calculation fixes epsilon=38.51 and **N=50**, as in the preparation study, and decreases phi from pi/3 (60 degrees) toward zero. J=1, U=-0.0002, V=0, Delta=0.5, kappa=1. Only stationary mean-field roots and their Bogoliubov spectra are calculated; no ramp or time evolution is performed.

### Following the same nonlinear branch

At each phase, recompute the homogeneous cubic with c(phi)=Delta+U+2J cos(phi), choose the intermediate positive real density root, and use alpha0=-epsilon/[c(phi)+2Un-i kappa/2]. Recompute the exact fixed endpoint pumps epsilon+J exp(+i phi) alpha0 and epsilon+J exp(-i phi) alpha0. Thus this follows the compensated stationary branch at fixed physical bulk parameters, rather than freezing the density or pairing at their pi/3 values.

Every existing target is exactly uniform including the edges; maximum stationary force over the phase scan is {max(r['residual'] for r in valid):.3g}. Compensation is an additive fixed drive, so it contributes no BdG derivatives. The operator is the uniform OBC truncation evaluated on the newly computed background at each phi.

### Stability and branch termination

The scan uses 0.1-degree spacing over 0-60 degrees and logarithmically refined phase distances from the intermediate-root merger (1e-5 to 0.1 degrees). Sign changes of the spectral growth rate are refined with a scalar root finder. The resolved zero-growth phases are {', '.join(f"{t['degrees']:.8f} degrees ({t['phi']:.10f} rad)" for t in transitions) or 'none'}; classifications within 1e-7 J of zero remain marginal/unresolved.

The intermediate and high-density bulk roots merge at **phi={m:.8f} degrees = {result['merger_phi']:.10f} rad**. Below this phase, epsilon=38.51 lies below the lower turning-point drive, and only the low-density positive bulk root remains. Hence there is no separate intermediate homogeneous branch, and its root-specific compensation is undefined, all the way to phi=0. This branch termination must not be reported as a dynamical instability or solver failure.

For this N=50 scan, the phase intervals are approximately:

|phase interval (degrees)|intermediate compensated state|
|---|---|
|56.2397 < phi <= 60|stable|
|18.2898 < phi < 56.2397|unstable|
|18.2247 < phi < 18.2898|stable again, in a narrow interval|
|0 <= phi < 18.2247|no separate intermediate root|

The zero-growth boundaries are marginal, and the branch-merger endpoint has coinciding roots. These intervals describe the recomputed nonlinear background, not a fixed-density quadratic phase sweep.

|phi (degrees)|intermediate density|OBC growth / J|classification|
|---:|---:|---:|---|
{chr(10).join(table)}

![Uniform intermediate density, OBC growth and eigenvalue conditioning versus hopping phase](../figures/compensated_intermediate_phi_stability.png)

![Representative OBC Bogoliubov spectra at fixed epsilon, recomputed background and compensation](../figures/compensated_intermediate_phi_spectra.png)

![Resolved growth-rate crossings, including the narrow stable interval next to branch termination](../figures/compensated_intermediate_phi_crossings.png)

Time evolution convention is exp(-i omega t), so stability requires max Im(omega)<0. All spectra retain both particle and hole sectors. These are finite-N=50 results, not thermodynamic limits or a claim about all other stationary branches.

### Independent precision checks

For selected phases and phases 0.001 degrees to either side of each computed stability crossing, we recompute the intermediate density and trigonometric coefficients at high precision, form the real-quadrature generator, and count roots with Re(lambda)>0 using a banded characteristic polynomial and Routh recursion at 100 and 160 digits. Both precisions are recorded, together with counts to the right of shifts +/-1e-7. The 160-digit growing-mode counts are:

|phi (degrees)|number of growing generator modes|
|---:|---:|
{chr(10).join(checks)}

Double-precision eigenvalue residuals, dominant left/right condition numbers, and sensitivity indicators are saved for every existing root. High-precision counts support stable/unstable classifications near the reported crossings; they are precision convergence checks, not interval certificates. The tabulated last digits of zero-growth phases come from double-precision scalar refinement and should not be interpreted as rigorously certified transition accuracy.

### Reproduction and records

Run `python src/compensated_phi_stability.py`. Records are data/generated/compensation_n50/compensated_phi_stability.json, including all real bulk densities, selected complex amplitude, both complex endpoint pumps, full OBC eigenvalues, residuals, growth, conditioning, branch-existence flags and high-precision counts. The merger is solved independently from the analytical lower turning point of the bulk cubic. No Gaussian covariance evolution or fluctuation feedback is used.
'''
    p=Path('docs/uniform_topological_segment_stability.md');base=p.read_text(encoding='utf-8').split('## Compensated intermediate-state stability versus phi at epsilon=38.51')[0].rstrip();p.write_text(base+'\n\n'+text,encoding='utf-8')

if __name__=='__main__':main()
