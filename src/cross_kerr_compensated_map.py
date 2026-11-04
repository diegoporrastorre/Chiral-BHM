"""N=50, U=V=-2e-4, Delta=0, fixed-drive compensated stability map."""
import compensated_stability_map as drawing
from compensated_stability_map import *
from decimal import Decimal,localcontext
from banded_stability_decimal import characteristic_banded,count_right
from scipy.optimize import linear_sum_assignment

V=U;DELTA=0.;beta=2*(U+V)
OUT=Path('data/generated/compensation_n50_u_equal_v_delta0');OUT.mkdir(parents=True,exist_ok=True)
CROSS_CUSP=float(np.arccos((np.sqrt(3)*gamma-DELTA-U)/2))
C=np.diag(np.ones(N-1),1)+np.diag(np.ones(N-1),-1);I=np.eye(N);degree=C.sum(axis=1)

def cross_bounds(phi):
    c=DELTA+U+2*np.cos(phi);rad=np.sqrt(np.maximum(c*c-3*gamma*gamma,0))
    return (-2*c+rad)/(3*beta),(-2*c-rad)/(3*beta)

def cross_drive(n,phi):return np.sqrt(n*((DELTA+U+2*np.cos(phi)+beta*n)**2+gamma*gamma))

def blocks(n,phi,size=N):
    adj=np.diag(np.ones(size-1),1)+np.diag(np.ones(size-1),-1)
    t=np.diag(np.full(size-1,np.exp(-1j*phi)),1)+np.diag(np.full(size-1,np.exp(1j*phi)),-1)
    a=t+V*n*adj+np.diag(DELTA+U+4*U*n+V*n*adj.sum(axis=1)-1j*gamma)
    b=2*U*n*np.eye(size)+V*n*adj
    return a,b

def cross_matrix(n,phi,size=N):
    a,b=blocks(n,phi,size);return np.block([[a,b],[-b.conj(),-a.conj()]])

def cross_growth(n,phi):return float(eigvals(cross_matrix(n,phi),check_finite=False).imag.max())

def stationary_model(n,phi):
    eps=float(cross_drive(n,phi));alpha=-eps/(DELTA+U+2*np.cos(phi)+beta*n-1j*gamma)
    m=Chain(N=N,drive=eps,phi=phi,U=U,V=V,delta=DELTA);m.pump=m.pump.astype(complex)
    m.pump[0]+=(np.exp(1j*phi)+V*n)*alpha
    m.pump[-1]+=(np.exp(-1j*phi)+V*n)*alpha
    return m,alpha

def decimal_check(degrees,epsilon,n_initial,precision):
    with localcontext() as ctx:
        ctx.prec=precision+20;D=Decimal;tol=D(10)**(-precision-10)
        def atan(x):
            total=term=x;k=1
            while abs(term)>tol:term=-term*x*x*(2*k-1)/(2*k+1);total+=term;k+=1
            return total
        pi=16*atan(D(1)/5)-4*atan(D(1)/239);phi=pi*D(str(degrees))/180
        cs=term=D(1);k=1
        while abs(term)>tol:term=-term*phi*phi/((2*k-1)*(2*k));cs+=term;k+=1
        sn=term=phi;k=1
        while abs(term)>tol:term=-term*phi*phi/((2*k)*(2*k+1));sn+=term;k+=1
        ud=D('-.0002');vd=ud;bd=2*(ud+vd);c=ud+2*cs;n=D(str(n_initial));e=D(str(epsilon))
        for _ in range(100):
            step=(n*((c+bd*n)**2+D('.25'))-e*e)/((c+bd*n)*(c+3*bd*n)+D('.25'));n-=step
            if abs(step)<tol:break
        ctx.prec=precision;cs=+cs;sn=+sn;n=+n
        gs=-2*ud*n;gc=-vd*n;ar=cs+vd*n
        a=[[D(0)]*(2*N) for _ in range(2*N)]
        for j in range(N):
            z=1 if j in [0,N-1] else 2
            d=ud+4*ud*n+z*vd*n
            a[2*j][2*j]=gs-D('.5');a[2*j+1][2*j+1]=-gs-D('.5')
            a[2*j][2*j+1]=d;a[2*j+1][2*j]=-d
            if j<N-1:
                for p,q,im in [(j,j+1,-sn),(j+1,j,sn)]:
                    a[2*p][2*q]=im+gc;a[2*p+1][2*q+1]=im-gc
                    a[2*p][2*q+1]=ar;a[2*p+1][2*q]=-ar
        coeff=characteristic_banded(a,3)
        return dict(precision=precision,N=N,degrees=float(degrees),epsilon=float(epsilon),density=float(n),
                    counts={str(s):count_right(coeff,D(str(s))) for s in [-1e-7,0,1e-7]})

def validate():
    # Direct model and independently built operator agree after a global pairing-phase rotation.
    rng=np.random.default_rng(20261104)
    for phi in [0.,.7,1.]:
        n=700.;model=Chain(N=8,V=V,U=U,delta=0,phi=phi);real=np.full(8,np.sqrt(n))
        np.testing.assert_allclose(cross_matrix(n,phi,8),model.bogoliubov(real),atol=1e-12)
        a=rng.normal(size=8)+1j*rng.normal(size=8);direction=rng.normal(size=8)+1j*rng.normal(size=8)
        aa,bb=model.blocks(a);h=1e-5
        err=max(abs((model.force(a+h*direction)-model.force(a-h*direction))/(2*h)-aa@direction-bb@direction.conj()))
        assert err<1e-8
        # Check the real-quadrature matrix used by the Decimal root count.
        A,B=blocks(n,phi,8);Bimag=-1j*B
        G=np.block([[(A+Bimag).imag,(A-Bimag).real],[-(A+Bimag).real,(A-Bimag).imag]])
        direct=-1j*eigvals(cross_matrix(n,phi,8));real_eig=eigvals(G)
        costs=abs(direct[:,None]-real_eig[None,:]);rows,cols=linear_sum_assignment(costs)
        assert max(costs[rows,cols])<1e-10

def main():
    validate()
    phase_edges=np.linspace(0,CROSS_CUSP,321);phi=(phase_edges[:-1]+phase_edges[1:])/2
    base=np.sin(np.linspace(0,np.pi/2,81))**2
    f=np.unique(np.r_[base,[1e-7,1e-6,1e-5,1e-4,.001,.003],1-np.array([1e-7,1e-6,1e-5,1e-4,.001,.003])]);f=f[(f>0)&(f<1)]
    f_edges=np.r_[0,(f[:-1]+f[1:])/2,1];a,b=cross_bounds(phi)
    density=a[:,None]+(b-a)[:,None]*f;epsilon=cross_drive(density,phi[:,None]);rates=np.full(density.shape,np.nan)
    path=OUT/'compensated_stability_map.npz'
    if path.exists():
        old=np.load(path)
        if np.array_equal(old['phi'],phi) and np.array_equal(old['fraction'],f):rates=old['growth'].copy()
    for i,p in enumerate(phi):
        if np.isfinite(rates[i]).all():continue
        rates[i]=[cross_growth(float(n),float(p)) for n in density[i]]
        if i%20==0 or i==len(phi)-1:
            np.savez_compressed(path,phi=phi,phase_edges=phase_edges,fraction=f,fraction_edges=f_edges,density=density,epsilon=epsilon,growth=rates)
            print('row',i,'of',len(phi),flush=True)
    curves=[]
    for i,p in enumerate(phi):
        for j in range(len(f)-1):
            if rates[i,j]*rates[i,j+1]<0:
                n=brentq(lambda n:cross_growth(n,p),density[i,j],density[i,j+1],xtol=1e-7)
                curves.append(dict(phi=float(p),epsilon=float(cross_drive(n,p)),density=float(n)))
    selected=[]
    for deg in [0,15,30,45,60]:
        p=np.deg2rad(deg);lo,hi=cross_bounds(p)
        for t in [.005,.5,.995]:
            n=float(lo+(hi-lo)*t);e=float(cross_drive(n,p));g=cross_growth(n,p)
            checks=[decimal_check(deg,e,n,precision) for precision in [100,160]]
            assert checks[0]['counts']==checks[1]['counts'];assert (g<0)==(checks[1]['counts']['0']==0)
            m,alpha=stationary_model(n,p);residual=float(max(abs(m.force(np.full(N,alpha)))))
            assert residual<1e-8
            selected.append(dict(degrees=deg,epsilon=e,density=n,growth=g,stationary_residual=residual,checks=checks,
                                 endpoint_pumps=[[float(m.pump[j].real),float(m.pump[j].imag)] for j in [0,N-1]]))
            print('precision',deg,e,g,checks[-1]['counts']['0'],flush=True)
    summary=dict(N=N,parameters=dict(J=1,U=U,V=V,delta=DELTA,kappa=1),cusp_phi=CROSS_CUSP,cusp_degrees=float(np.rad2deg(CROSS_CUSP)),
                 samples=int(rates.size),stable_samples=int((rates<-1e-7).sum()),unstable_samples=int((rates>1e-7).sum()),marginal_samples=int((abs(rates)<=1e-7).sum()),
                 boundary_points=curves,precision_checks=selected)
    (OUT/'compensated_stability_map.json').write_text(json.dumps(summary,indent=2))
    drawing.N=N;drawing.DEST=OUT;drawing.CUSP=CROSS_CUSP;drawing.bounds=cross_bounds;drawing.drive=cross_drive
    drawing.FIGURE_PREFIX='compensated_stability_map_u_equal_v_delta0'
    drawing.PLOT_TITLE='Compensated homogeneous intermediate state: OBC stability, N=50\nU=V=-0.0002, Delta=0; J=kappa=1'
    drawing.PLOT_RANGES=[('full',(0,85,0,np.pi)),('zoom',(10,45,0,np.deg2rad(70)))]
    drawing.plot(phi,phase_edges,f_edges,rates,curves);report(summary)

def report(summary):
    cusp_e=np.sqrt((-2*np.sqrt(3)*gamma/(3*beta))*(gamma*gamma/3+gamma*gamma));lo,hi=cross_bounds(0)
    text=f'''## Cross-Kerr stability map: N=50, U=V=-0.0002, Delta=0

This repeats the compensation-stabilizability calculation with **N=50, U=V=-2e-4 and Delta=0**, retaining J=1 and kappa=1. Axes remain epsilon horizontally and phi vertically, over 0<=phi<=pi. No ramp or time evolution is calculated. Green means an exact homogeneous compensated stationary solution whose actual cross-Kerr OBC BdG operator is stable; orange means unstable; grey means no distinct intermediate positive-density root; yellow means marginal/unresolved growth within 1e-7 J.

![Full compensation-stability map with equal onsite and nearest-neighbour interactions, Delta=0](../figures/compensated_stability_map_u_equal_v_delta0_full.png)

![Zoom on the complete intermediate-root lobe for equal U and V, Delta=0](../figures/compensated_stability_map_u_equal_v_delta0_zoom.png)

### Bulk root and the necessary change in compensation

Each nearest-neighbour cross-Kerr bond is counted once, H_cross=V sum_j n_j n_(j+1). The normally ordered mean-field onsite term retains the +U shift from U n_j^2. A uniform bulk root obeys

`epsilon^2=n[(Delta+U+2J cos(phi)+2(U+V)n)^2+(kappa/2)^2]`.

Here the nonlinear coefficient in this cubic is 2(U+V)=-0.0008; c(phi)=-0.0002+2 cos(phi). The middle of three positive density roots defines the intermediate solution. Its complex alpha0 is -epsilon/[c(phi)+2(U+V)n-i kappa/2].

An open end lacks both a hopping contribution and one cross-Kerr mean-field contribution. Keeping **Delta=0 at every site**, the required fixed additive endpoint pumps are

`pump_1=epsilon+[J exp(+i phi)+V n] alpha0`,

`pump_N=epsilon+[J exp(-i phi)+V n] alpha0`.

All other pumps remain epsilon. Omitting the V n alpha0 term would fail to make the state homogeneous and stationary. We compensate through fixed coherent drives, following the earlier protocol; no endpoint detuning or feedback is introduced. This prescription is root-specific, not optimized.

### The cross-Kerr BdG operator is not an onsite-Kerr substitution

For a uniform alpha0 and open-chain adjacency C, let z_j=sum_k C_jk (1 at each end, 2 in the interior). The exact linearization is

`A_jk=T_jk+delta_jk[Delta+U+4Un+z_j Vn-i kappa/2]+Vn C_jk`,

`B_jk=2U alpha0^2 delta_jk+V alpha0^2 C_jk`,

`L=[[A,B],[-B*,-A*]]`.

Thus cross-Kerr coupling produces both normal and anomalous intersite terms. In addition, the onsite coefficient at an end differs from the bulk by -Vn because that site has only one neighbour. **Fixed pump compensation cancels the stationary mismatch but does not restore missing fluctuation couplings or their onsite derivatives.** This map uses those actual degree-one edge coefficients. Replacing U by U+V in the earlier onsite-only BdG matrix, or silently using degree two at the ends, would describe a different operator and give a different map. Compensating through endpoint detunings would also be a different boundary prescription.

Modes evolve as exp(-i omega t); max Im(omega)<0 gives local stability. A global Nambu rotation makes the common pairing phase convenient without changing any eigenvalue. Because the compensation is fixed and additive, it contributes no BdG derivatives.

### Existence region, resolution and checks

Positive-density bistability requires c(phi)>sqrt(3) kappa/2, restricting the intermediate-root lobe to phi<{summary['cusp_degrees']:.6f} degrees. The cusp drive is epsilon={cusp_e:.8f}; at phi=0 the root-existence endpoints are epsilon={float(cross_drive(hi,0)):.8f} and {float(cross_drive(lo,0)):.8f}. Existence curves exclude the merger points as distinct intermediate roots. The full plot uses the earlier 0-85 drive scale for comparison, while the zoom covers 10-45.

The curved sampling has 320 phase rows and 91 intermediate-density fractions, clustered near both root mergers: {summary['samples']} spectra, of which {summary['stable_samples']} are stable, {summary['unstable_samples']} unstable and {summary['marginal_samples']} marginal. These counts are not geometric areas. Growth-rate zeros along each row are refined; the shading interpolates the spectra with those zero-growth points and masks outside the analytic existence lobe. The dashed vertical line marks epsilon=38.51 for comparison with earlier results. Near fine boundaries the plot is a finite-resolution numerical characterization.

The independent matrix construction is checked against the reusable cross-Kerr mean-field implementation, and its analytic Jacobian against finite differences. Representative phases 0,15,30,45,60 degrees at three density fractions are checked with the actual compensated complex stationary force (all below 1e-8). At those points independent real-quadrature characteristic-polynomial/Routh counts at 100 and 160 digits agree with the double-precision stable/unstable labels. These counts include onsite degree-one edge terms and both intersite pairing and normal coupling. Selected precision checks are not rigorous certificates for every plotted point.

Green regions establish locally stable compensated mean-field states only, not successful preparation from vacuum or robustness to compensation errors. The earlier V=0 maps remain in the report for comparison, but both Delta and V have changed here, so differences cannot be attributed to V alone. No Gaussian covariance evolution or fluctuation feedback enters this calculation.

Reproduce with `python src/cross_kerr_compensated_map.py`. Grid arrays, densities, drives, phases and growth rates are stored in data/generated/compensation_n50_u_equal_v_delta0/compensated_stability_map.npz. The companion JSON records parameters, refined zero-growth curves, endpoint pumps, stationary residuals and 100/160-digit checks. Full and zoom PNG/PDF figures use prefix compensated_stability_map_u_equal_v_delta0 under figures/.
'''
    p=Path('docs/uniform_topological_segment_stability.md');base=p.read_text(encoding='utf-8').split('## Cross-Kerr stability map: N=50, U=V=-0.0002, Delta=0')[0].rstrip();p.write_text(base+'\n\n'+text,encoding='utf-8')

if __name__=='__main__':main()
