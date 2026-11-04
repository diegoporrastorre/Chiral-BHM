"""Add endpoint compensation experiment to the existing report."""
from pathlib import Path
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]

def main():
    r=json.loads((ROOT/'data/generated/homogeneous_obc_n100/endpoint_compensation.json').read_text())
    table=[]
    for b in r['branches']:
        first=b['continuation'][0];last=b['continuation'][-1]
        table.append(f"|{b['label']}|{b['uniform_density']:.3f}|{first['growth']:.6f}|{first['smin']:.4g}|{b['initial_static_susceptibility']['max_dalpha_ds']:.4g}|{b['continuation_end_s']:.9f}|")
    mid=r['branches'][1];small=[]
    for target in [1e-9,1e-8,1e-7,1e-6,.01]:
        x=min(mid['continuation'],key=lambda x:abs(1-x['s']-target))
        small.append(f"|{1-x['s']:.4g}|{np.mean(x['density']):.3f}|{x['growth']:.6f}|")
    section=f'''## Endpoint-compensation experiment at epsilon=38.51

This experiment separates stationary-background selection from BdG stability. It supports the interpretation that a stable topological background can be strongly deformed by boundary forcing; it does not establish a general front-selection theory.

### Construction

The three homogeneous bulk roots satisfy

`epsilon^2 = n [(Delta+U+2J cos(phi)+2Un)^2+(kappa/2)^2]`,

with `alpha0=-epsilon/(Delta+U+2J cos(phi)+2Un-i kappa/2)`.

For each root independently, set the additive complex endpoint drives to

`pump_1(s)=epsilon+s J exp(+i phi) alpha0`,

`pump_N(s)=epsilon+s J exp(-i phi) alpha0`.

Interior drives are unchanged. The complex drive Hamiltonian is pump_j a_j^dagger + pump_j* a_j. These fixed external terms replace the mean field of the missing left/right neighbours at s=1; they do not couple to fluctuations. At s=1 the uniform profile is therefore an exact OBC stationary solution. At s=0 the model is exactly the original uniformly driven open chain. The different bulk roots define different endpoint-drive paths for 0<s<=1, although all paths end at the same model at s=0.

The code verifies initial stationary residuals below 1e-10 and exact equality of compensated/uncompensated BdG matrices evaluated at the same uniform background. Along the path the background changes, so its BdG matrix generally changes too.

### Results: stability versus boundary sensitivity

|starting root|uniform density|growth at s=1 / J|smallest singular value at s=1 / J|max abs(d alpha_j/ds) at s=1|last locally continued s|
|---|---:|---:|---:|---:|---:|
{chr(10).join(table)}

All three starting uniform OBC backgrounds are linearly stable. Nevertheless, the intermediate topological root has a dramatically larger response to the particular missing-neighbour forcing. This is not inferred just from its small singular value: we explicitly solve

`L z' = -(q,-q*)`,

where q_1=J exp(+i phi) alpha0 and q_N=J exp(-i phi) alpha0, other entries vanish, and the first N components of z' give d alpha/ds. The intermediate result max abs(d alpha/ds)=1.15e9 means that tiny changes of compensation can cause substantial deformation. It is a local derivative, not a finite-amplitude prediction outside linear response.

The intermediate branch can be continued beyond the uniform point if compensation steps begin at 1e-11. The original minimum step of about 1e-7 was inadequate; its failure did not indicate absence of a stationary branch. The refined branch shows the following actual nonlinear changes:

|removed compensation 1-s|mean density|growth / J|
|---:|---:|---:|
{chr(10).join(small)}

![Intermediate-root sensitivity resolved near full compensation](../figures/endpoint_compensation_middle_refinement.png)

Thus **substantial inhomogeneity develops while the intermediate branch remains stable at the displayed points**. A dynamical instability is not necessary to obtain that deformation. The full OBC operator has large static susceptibility despite decaying eigenmodes.

### Continuation and dynamics are separate observations

![Growth along stationary continuation and selected dynamically reached density profiles](../figures/endpoint_compensation_paths.png)

The high-root stationary branch is followed continuously from s=1 to s=0 and stays stable: its largest computed growth along the path is {max(x['growth'] for x in r['branches'][2]['continuation']):.6f} J. It ends at the original inhomogeneous profile. This is the clearest example of a stable homogeneous background deforming into the observed nonuniform state without a resolved stability loss.

Local continuation from the low root stops near s={r['branches'][0]['continuation_end_s']:.9f}, with the final computed growth {r['branches'][0]['continuation'][-1]['growth']:.3g} J. The refined intermediate path stops near s={mid['continuation_end_s']:.9f}, with growth {mid['continuation'][-1]['growth']:.3g} J. These near-marginal endpoints are consistent with branch turning or bifurcation, but the stopping rule is numerical: root acceptance fails or the permitted profile step becomes too small. Pseudo-arclength continuation is needed to establish the branch geometry. No fold is claimed as proven. Growth magnitudes below 1e-7 J are treated as marginal/unresolved; the intermediate path has no resolved positive growth above that threshold.

Separately, all three initial conditions were evolved under a staircase protocol s=1,0.98,...,0. Each step runs to t=1200 with DOP853, rtol=2e-10, atol=2e-11, starting from the preceding endpoint. All 153 steps finish close to root-polished stationary solutions (max amplitude distance <1e-3 and max force <1e-3). The reported profiles use these reached stationary roots; the integration endpoint, final-50-time density variation, and root distance are recorded separately. This protocol has finite drive steps and does not locate bifurcations precisely. In particular, its first intermediate-root step skips the extremely narrow near-s=1 deformation scale.

![Density evolution under the staircase compensation protocol](../figures/endpoint_compensation_density_maps.png)

At s=0 all three protocols recover the original epsilon=38.51 profile, with maximum site-density differences {', '.join(f"{b['final_density_difference_from_original']:.3g}" for b in r['branches'])} for low/middle/high starts respectively. This is evidence for a common selected endpoint under these protocols, not proof of uniqueness or absence of other basins.

### What this establishes

The compensation experiment weakens the explanation that spatial variation is required only as a stabilizing taper. It supplies a concrete alternative: open boundaries impose a finite stationary mismatch, and a linearly stable topological background can respond extremely strongly to that mismatch. The nonlinear stationary equations then determine a deformed profile, potentially connecting different bulk density branches. The high-root continuation demonstrates stable deformation all the way to the target model; the intermediate-root response demonstrates the extreme sensitivity most directly.

The results do not yet establish that topology is the unique cause of the sensitivity, a universal relation between topology and profile selection, or an analytic law for interface position. Comparisons with phi=0 or reversed phi, system-size scaling, and pseudo-arclength continuation would address those remaining questions. The susceptibility calculation is specific to the physical endpoint-forcing vector, which is more informative here than the norm of L inverse alone.

### Reproduce and numerical limits

```powershell
python src/endpoint_compensation.py
python src/refine_endpoint_compensation.py
python src/report_endpoint_compensation.py
```

All parameters remain N=100, J=1, U=-0.0002, V=0, Delta=0.5, kappa=1, phi=pi/3. The bulk drive is fixed at epsilon=38.51. Stationary root tolerance is 1e-10, acceptance max|f|<1e-8. Local continuation uses secant prediction, maximum compensation step 0.005, step halving after failure or excessive amplitude change; the initial runs reject changes above 2 and stop below step 1e-7, while the intermediate refinement rejects changes above 1 and stops below 1e-13. This is parameter continuation, not pseudo-arclength continuation.

All complex amplitudes, densities, stationary growth, residuals, singular diagnostics, failed steps, and time-evolution convergence diagnostics are saved in data/generated/homogeneous_obc_n100/endpoint_compensation.json. The near-uniform intermediate operator is highly nonnormal; calculations here use double precision. Its previously checked stability margin is substantial, but tiny near-marginal rates and late continuation endpoints are not certified high-precision thresholds. The large measured susceptibility and resolved finite deformations support sensitivity without relying on the last digits of those thresholds.
'''
    p=ROOT/'docs/uniform_topological_segment_stability.md'
    base=p.read_text(encoding='utf-8').split('## Endpoint-compensation experiment at epsilon=38.51')[0].rstrip()
    p.write_text(base+'\n\n'+section,encoding='utf-8')

if __name__=='__main__':main()
