"""Append the N=50 compensation and vacuum-start results."""
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]

def main():
    dest=ROOT/'data/generated/compensation_n50'
    result=json.loads((dest/'compensation.json').read_text());scan=json.loads((dest/'scan_35_45.json').read_text())
    difference=max(max(abs(np.array(a['density'])-np.array(b['density']))) for a,b in zip(scan[:1001],reversed(scan[1001:])))
    table=[];pumps=[]
    for b in result['branches']:
        z=b['vacuum'];zz=b['vacuum_refined']
        table.append(f"|{b['label']}|{b['density_target']:.3f}|{b['target']['growth']:.6f}|{'yes' if z['target_reached'] else 'NO'}|{'yes' if zz['target_reached'] else 'NO'}|{'yes' if b['perturbed_target']['target_reached'] else 'NO'}|")
        f=lambda pair:f'{pair[0]:.6f} {pair[1]:+.6f}i'
        pumps.append(f"|{b['label']}|{f(b['endpoint_pumps'][0])}|{f(b['endpoint_pumps'][1])}|")
    b=result['branches'][1];z=b['vacuum'];zz=b['vacuum_refined']
    profile_diff=max(abs(np.array(z['density'])-np.array(zz['density'])))
    block=f'''## N=50: are compensated uniform states reachable from vacuum?

**For a sudden application of the fixed compensated drives, the low and high targets are reached, but the intermediate target is not.** All three targets remain exact, linearly stable OBC stationary solutions. Thus making a chosen root stationary and locally stable does not ensure that vacuum lies in its basin of attraction.

### Repeated baseline calculation

We repeat the homogeneous-drive OBC scan for N=50, epsilon=35 to 45 with spacing 0.01 in both directions, retaining J=1, U=-0.0002, V=0, Delta=0.5, kappa=1, phi=pi/3 and the U n^2 convention. Each direction begins from vacuum at its endpoint; subsequent points use stationary continuation and, when needed, evolution to obtain a new seed. This is not an independent vacuum-start calculation at every scan point.

All 2002 scan points converge with max residual {max(x['residual'] for x in scan):.3g}; all computed full-chain growth rates are negative (largest {max(x['growth'] for x in scan):.6f} J). Maximum up/down site-density difference is {difference:.3g}. The figures show the N=50 density profiles and local winding using the earlier zero-frequency convention. Local-bulk interpretation retains the slow-variation limitations stated above.

![N=50 fine baseline density and local winding maps](../figures/compensation_n50_scan.png)

### Fixed compensated systems at epsilon=38.51

For each of the same three bulk roots, compensation is recomputed from that root's complex alpha0. Only the two endpoints are modified. The complete additive endpoint drives (including the real bulk drive 38.51) are:

|target|pump at site 1 / J|pump at site 50 / J|
|---|---|---|
{chr(10).join(pumps)}

All interior additive drives remain 38.51. These are three different fixed Hamiltonians. At s=1 their respective uniform profiles satisfy the stationary equations with residual below 1e-10. Independent banded characteristic-polynomial/Routh counts at 100 and 160 digits find zero growing modes for each target. No covariance dynamics or fluctuation feedback is included.

### Direct evolution from zero photons

For each fixed compensated Hamiltonian separately, set **alpha_j(0)=0 at every site** and switch on the complete constant drives at t=0. There is no compensation ramp, target seeding, feedback, or intermediate root solve during this evolution. Integrate with DOP853 to t=1200 (rtol=2e-10, atol=2e-11). Repeat independently from vacuum to t=2400 with rtol=2e-12, atol=2e-13.

After integration only, root polish the endpoint to identify the approached stationary solution. A state is classified as stationary only if the polished residual is below 1e-8, the raw endpoint force is below 1e-3, and the maximum complex amplitude distance between the integration endpoint and polished root is below 1e-3. The target is reached only if, additionally, max_j|alpha_j-alpha0|<1e-3. Therefore a solver finding the desired root cannot substitute for actual dynamical convergence to it.

|target|uniform density|target BdG growth / J|vacuum reaches target, t=1200|vacuum reaches target, tighter t=2400|small perturbation returns to target|
|---|---:|---:|---|---|---|
{chr(10).join(table)}

The local-attraction check starts at alpha0 plus a complex random perturbation of scale 1e-4 per real/imaginary component, seed 20261104, and evolves under the same fixed drives to t=1200. All three return to the target. This tests one small perturbation, not an entire basin.

![N=50 direct vacuum dynamics and final profiles, compared with each intended uniform target](../figures/compensation_n50_vacuum.png)

For the intermediate compensation, vacuum instead approaches a stable **inhomogeneous** stationary state: mean density {np.mean(z['density']):.3f}, site-density range {min(z['density']):.3f}-{max(z['density']):.3f}, growth {z['growth']:.6f} J. Its maximum complex amplitude difference from the uniform target is {z['target_amplitude_error']:.3f}. The tighter, longer vacuum run reaches the same density profile to {profile_diff:.3g} maximum site-density difference. Its final raw force is {zz['end_force']:.3g}, and the maximum temporal density standard deviation in the last 50 time units is {zz['tail_density_std']:.3g}. Thus this failure to reach the intermediate target is not just incomplete relaxation at t=1200 or a loose integration tolerance.

The intermediate compensated system consequently supports **at least two distinct stable stationary states** in this mean-field calculation: the uniform target recovered from a nearby initial condition, and the inhomogeneous attractor reached from vacuum. Compensation stabilizes the target locally, but does not force vacuum to select it. We have not established whether a suitably designed preparation ramp could reach the intermediate target from vacuum; the negative result applies to the specified constant-drive switch-on protocol.

### Repeat removal of compensation at N=50

We also repeat local stationary continuation from s=1 toward zero for each root, plus the earlier s=1,0.98,...,0 staircase dynamics initialized at its uniform target. This preparation is deliberately distinguished from the direct-vacuum test above.

![N=50 compensation continuation growth and staircase density profiles](../figures/compensation_n50_continuation.png)

Local continuation reaches s={result['branches'][0]['continuation_end_s']:.9f} from the low root, s={result['branches'][1]['continuation_end_s']:.9f} from the intermediate root, and s={result['branches'][2]['continuation_end_s']:.9f} from the high root. The high branch reaches zero compensation. The other branches terminate numerically near marginal growth; these are not certified folds without pseudo-arclength continuation. Stationary continuation does not use a time-evolution fallback to jump to another branch.

### Reproduction, records, and limits

```powershell
python src/compensation_n50_vacuum.py
python src/report_compensation_n50.py
```

Driver: src/compensation_n50_vacuum.py. Generated data: data/generated/compensation_n50/scan_35_45.json and compensation.json. Records include target amplitudes and endpoint drives, BdG and high-precision counts, full vacuum trajectories, tighter-run trajectories, perturbed-target trajectories, raw endpoint force, root distance, final temporal variation, local continuation and failed attempts, and staircase profiles. Session interpreter is C:/Users/diego/miniconda3/envs/topo/python.exe.

Root tolerance is 1e-10, residual acceptance 1e-8. Stability tolerance is 1e-7 J. Local compensation continuation uses secant prediction, maximum step 0.005, step halving after failed roots or maximum amplitude changes above 1; minimum step is 1e-12, initial intermediate step 1e-7. The baseline scan uses the same mean-field-only evolution fallback as before. The sharper near-marginal continuation endpoints are numerical observations, not high-precision phase boundaries. This experiment concerns N=50 at epsilon=38.51; it does not assert vacuum reachability for every drive or chain length.
'''
    p=ROOT/'docs/uniform_topological_segment_stability.md';text=p.read_text(encoding='utf-8').split('## N=50: are compensated uniform states reachable from vacuum?')[0].rstrip()
    p.write_text(text+'\n\n'+block,encoding='utf-8')

if __name__=='__main__':main()
