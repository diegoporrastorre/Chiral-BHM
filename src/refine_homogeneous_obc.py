"""Refine the extended topological domain and write the scan report."""
from scan_homogeneous_obc import *
from bdg_topology import winding

def runs(mask):
    ids=np.where(mask)[0]+1
    return ', '.join(str(g[0]) if len(g)==1 else f'{g[0]}-{g[-1]}' for g in np.split(ids,np.where(np.diff(ids)>1)[0]+1) if len(g)) or 'none'

def main():
    rr=json.loads((OUT/'scan.json').read_text()); up=rr[:161]; down=rr[161:]
    seed=next(r for r in up if r['epsilon']==38.25)
    a=np.array(seed['alpha_real'])+1j*np.array(seed['alpha_imag'])
    refined=[]
    for eps in np.arange(38.25,38.801,.01):
        m=Chain(N=100,drive=float(eps));r=m.stationary(a,1e-8)
        if not r['converged']:
            a,_,_=evolve(m,a,duration=1200,rtol=2e-10,atol=2e-11)
            r=m.stationary(a,1e-8)
        assert r['converged'],eps
        a=r['alpha'];rec=diagnose(m,a);rec['epsilon']=float(eps)
        refined.append(rec)
    (OUT/'refinement.json').write_text(json.dumps(refined,indent=2))
    # Check analytic winding against momentum-space determinant at every site
    # of the most extended domain. Sign follows Fourier exp(+ikj).
    best=max(refined,key=lambda r:sum(np.array(r['nu'])!=0))
    a=np.array(best['alpha_real'])+1j*np.array(best['alpha_imag']);n=abs(a)**2
    for j in range(100):
        w=winding(1,np.pi/3,.5-2e-4-8e-4*n[j],-4e-4*a[j]**2,.5,nk=16384)
        assert w['winding']==best['nu'][j],(j,w)
    checks=json.loads((OUT/'checks.json').read_text())
    difference=max(max(abs(np.array(u['density'])-np.array(d['density']))) for u,d in zip(up,reversed(down)))
    vacuum_diff=max(max(abs(np.array(c['density'])-np.array(next(r for r in up if r['epsilon']==c['epsilon'])['density']))) for c in checks)
    # Save compact CSV profiles independently of JSON format.
    np.savetxt(OUT/'density_up.csv',np.column_stack(([r['epsilon'] for r in up],np.array([r['density'] for r in up]))),delimiter=',',header='epsilon,'+','.join(f'n{j}' for j in range(1,101)),comments='')
    rows=[]
    for eps in [20,30,35,38,38.5,38.75,40,42,45,48,50,60]:
        r=next(x for x in up if x['epsilon']==eps);n=np.array(r['density'])
        rows.append(f"| {eps:g} | {n.mean():.2f} | {n.min():.2f}-{n.max():.2f} | {r['growth']:.5f} | {runs(np.array(r['nu'])!=0)} | {r['smin']:.3g} |")
    q=np.array(best['nu'])!=0; eta=np.array(best['variation'])
    fig,ax=plt.subplots(3,1,figsize=(10,8),layout='constrained')
    ax[0].plot(range(1,101),best['density']);ax[0].set(ylabel='density',title=f"epsilon={best['epsilon']:.2f}; N=100, homogeneous OBC")
    ax[1].step(range(1,101),best['nu'],where='mid');ax[1].set(ylabel='local winding',yticks=[-1,0])
    ax[2].semilogy(range(1,101),eta);ax[2].axhline(.05,ls='--',color='k');ax[2].set(ylabel='bond variation / J',xlabel='site j')
    fig.savefig('figures/homogeneous_obc_n100_topological_profile.png',dpi=180)
    text=f'''# Homogeneous driven OBC chain: N=100, Delta=0.5

Calculated 2026-11-04. This is a mean-field stationary-state calculation followed by Bogoliubov analysis, with no Gaussian covariance evolution or feedback. These are new numerical observations, not a reproduction of the paper's Gaussian phase diagram.

## Parameters and conventions

J=1, U=-0.0002, kappa=1, phi=pi/3, Delta=0.5, V=0, N=100; epsilon runs from 20 to 60. All rates are in units J and time in 1/J. Drive and detuning are uniform on every site, including both edges. Sites are 1 through 100; only bonds (j,j+1) are present. In the real-drive gauge T[j+1,j]=J exp(+i phi). Homogeneous parameters do not require a homogeneous stationary density under OBC.

The Hamiltonian onsite convention is U n_j^2, so normal ordering retains the +U shift:

`i d alpha_j/dt = epsilon + (Delta+U-i kappa/2+2U|alpha_j|^2) alpha_j + J exp(-i phi) alpha_(j+1) + J exp(+i phi) alpha_(j-1)`.

Missing neighbours are omitted. Density means n_j=|alpha_j|^2. We solve the complex stationary force f=0. Fluctuations obey i dz/dt=L z, L=[[A,B],[-B*, -A*]], A=T+diag(Delta+U+4Un-i kappa/2), B=diag(2U alpha^2). Modes evolve as exp(-i omega t); stability therefore requires max Im(omega)<0 (classification tolerance 1e-7).

## Results

All 161 sampled drives (spacing 0.25) admit a converged, linearly stable stationary solution in both continuation directions. Maximum residual max_j|f_j| is {max(r['residual'] for r in rr):.3g}. The least negative computed growth rate is {max(r['growth'] for r in rr):.6f}. Up/down profiles agree to {difference:.3g} in maximum absolute density difference; distinct stable profiles occur at epsilon=47.25 and 47.50, providing evidence of a small hysteretic or multistable interval. The other sampled drives agree within 0.01 in absolute density. This does not enumerate all stationary branches.

Low drives have a low-density interior with appreciable boundary distortions. Near epsilon=38.5 an extended, slowly varying intermediate-density domain develops. Beyond this narrow crossover the interior is high density, with a low-density segment and an interface towards the left end. Increasing drive moves the interface leftwards. Even at epsilon=60 the end-site density is strongly different from the bulk. A description in terms of uniform density phases alone misses these spatial structures.

| epsilon | mean density | min-max density | max Im omega | nonzero frozen local winding: sites | smallest singular value |
|---:|---:|---:|---:|---|---:|
{chr(10).join(rows)}

![Density and local winding maps](../figures/homogeneous_obc_n100_scan.png)

![Selected profiles and stability](../figures/homogeneous_obc_n100_profiles.png)

## Topology: definition and interpretation

We explicitly checked the main paper and Section VI, Eqs. (15)-(22), of the [accessible arXiv v1 supplemental text](https://arxiv.org/html/2411.08965v1#S6). This is not a claim to have inspected a separate published supplement. The published main source is [Rassaert et al., PRL 135, 203603 (2025)](https://doi.org/10.1103/w32k-x3l5), stored in literature/papers/RassaertPRL2025.pdf.

At reference frequency omega_0=0, Q=-L, and H_ext=[[0,Q],[Q^dagger,0]]. Its chiral symmetry operator is diag(I,-I); its energies are plus/minus the singular values of Q. For the frozen homogeneous symbol at site j define d_j=Delta+U+4Un_j, b_j=2U alpha_j^2. With Fourier convention exp(+ikj),

`L_j(k)=(2J sin(phi) sin(k)-i kappa/2) I + [[d_j+2J cos(phi) cos(k), b_j],[-b_j*, -d_j-2J cos(phi) cos(k)]]`.

The local diagnostic is nu_j=(1/2pi) Delta_arg det[-L_j(k)] over increasing k from -pi to pi. Let X0=|b_j|^2-(kappa/2)^2-(d_j+2J cos(phi))^2 and Xpi=|b_j|^2-(kappa/2)^2-(d_j-2J cos(phi))^2. Nonzero winding occurs exactly when X0 Xpi<0; its sign is sign(d_j) for this gauge. Thus nontrivial sites here have nu=-1; reversing orientation reverses the sign. At either zero the point gap closes and the winding is undefined. The phase of b does not change the frozen invariant, but its spatial variation matters for using it locally.

The main extended topological candidate is near epsilon=38.5, rather than every isolated coloured pixel in the local map. On the refined 0.01 grid (38.25-38.80), the largest domain occurs at epsilon={best['epsilon']:.2f}, with nonzero winding on sites {runs(q)}. Sites {runs(q & (eta<.05))} also satisfy the descriptive slow-variation check max(|d_(j+1)-d_j|,|b_(j+1)-b_j|)<0.05J on adjacent bonds. This threshold is a diagnostic, not a rigorous error bound. Near a local gap closing even small variation can invalidate a frozen-bulk interpretation. At the selected point all 100 analytic winding labels agree with a 16384-point numerical momentum winding.

![Extended local topological domain and spatial variation](../figures/homogeneous_obc_n100_topological_profile.png)

At this point s_min={best['smin']:.6g}, s_next={best['snext']:.6g}; the right/left smallest-singular-vector centres are {best['singular_centers'][0]:.2f} and {best['singular_centers'][1]:.2f}. This provides a finite-chain SVD response diagnostic for the inhomogeneous operator, independently of the local approximation. It does not by itself establish protected modes or asymptotic exponential scaling. Most one- or two-site nonzero regions at other drives occur in rapidly varying boundary layers and cannot be identified confidently as local bulk topological phases. Topology and dynamical stability remain separate diagnostics.

## Numerical method, checks, and limits

Root solving uses the analytic real Jacobian, scipy.optimize.root with tol=1e-10, and acceptance max|f|<1e-8. Each sweep begins from vacuum at its endpoint; subsequent points use the previous profile. If a root fails or is unstable, mean-field DOP853 evolution over t=400 (rtol=2e-8, atol=2e-9) supplies another seed. A failed solve is never classified as absence of a stationary state. The refined scan continues from epsilon=38.25 using 0.01 steps, with a t=1200 tighter-tolerance evolution fallback when a root fails.

Independent vacuum starts at epsilon=20,30,35,38,40,42,45,50,60 were evolved to t=1200 with rtol=2e-10, atol=2e-11 and root polished. They agree with the up-sweep densities to {vacuum_diff:.3g}; polished residuals are below {max(c['residual'] for c in checks):.3g}. The largest temporal density standard deviation over the final 50 time units is {max(c['temporal_std'] for c in checks):.3g}, consistent with numerical integration error. The existing validation driver independently verifies the Jacobian against finite differences and the PBC spectrum against its Bloch form.

Spectra and SVD use double precision. No thermodynamic phase boundary, size scaling, proof of uniqueness, or exhaustive search for unstable/disconnected branches is claimed. The narrow topological crossover is resolved at 0.01 only within the stated refinement window; the wider scan has 0.25 resolution. N is held at the requested 100, so finite-size dependence is not established. A nonzero frozen invariant in a short, rapidly varying region is reported as a candidate diagnostic, not a proven topological domain.

## Reproducibility and files

Run from the project root with numpy, scipy, matplotlib:

```powershell
python src/scan_homogeneous_obc.py
python src/refine_homogeneous_obc.py
python src/validate_mean_field.py
```

This session used C:/Users/diego/miniconda3/envs/topo/python.exe. Reusable drivers are in src/. Full parameters are fixed explicitly in the code and above. Full complex amplitudes, densities, residuals, growth rates, local winding, variation diagnostics, and singular data are in data/generated/homogeneous_obc_n100/scan.json and refinement.json; independent-start checks are in checks.json. density_up.csv contains every up-sweep density profile. Generated datasets are excluded from git; rerun the drivers to regenerate them.
'''
    Path('docs/homogeneous_obc_n100_results.md').write_text(text,encoding='utf-8')
    print('best',best['epsilon'],'sites',runs(q),'smin',best['smin'],'up/down',difference,'vacuum',vacuum_diff)

if __name__=='__main__':main()
