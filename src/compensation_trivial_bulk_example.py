"""Document an existing N=50 preparation result with a trivial uncompensated bulk."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from bdg_topology import winding


def main():
    source=Path('data/generated/compensation_n50/preparation_30_35/epsilon_32.000000.json')
    r=json.loads(source.read_text());b=r['baseline'];t=r['target']
    # Check analytic frozen labels independently with a finer momentum grid.
    for n,nu in zip(b['density'],b['nu']):
        assert winding(1,np.pi/3,.4998-.0008*n,.0004*n,.5,nk=16384)['winding']==nu
    n=r['density_target'];w=winding(1,np.pi/3,.4998-.0008*n,.0004*n,.5,nk=16384)
    assert w['winding']==-1 and all(run['hold']['target_reached'] for run in r['runs'])
    fig,ax=plt.subplots(figsize=(8,4),layout='constrained')
    ax.plot(range(1,51),b['density'],'o-',ms=3,label='uncompensated: vacuum attractor')
    ax.plot(range(1,51),t['density'],label='compensated: prepared intermediate branch')
    ax.set(xlabel='site j',ylabel='density',title='N=50, epsilon=32, phi=pi/3');ax.legend()
    fig.savefig('figures/compensation_trivial_bulk_example.png',dpi=180)
    text='''## Example: compensation creates an accessible topological bulk at epsilon=32

At N=50, J=kappa=1, U=-0.0002, V=0, Delta=0.5, phi=pi/3 and epsilon=32, the existing independent vacuum-preparation calculation gives the following comparison. Stability uses i*dot(z)=L*z and the spectral reference for winding/SVD is omega=0.

| Quantity | Uncompensated vacuum attractor | Compensated intermediate state |
|:---|---:|---:|
| Density at site 25 | 535.695 | 2988.597 |
| Density range, sites 10–40 | 468.119–616.593 | 2988.597 (uniform) |
| Frozen bulk winding at site 25 | 0 | -1 |
| max Im(omega)/J | -0.261254 | -0.102268 |
| 1/s0 in 1/J | 5.32683 | 37963.3 |
| Stationary max force/J | 1.39e-14 | 4.94e-14 |

![Uncompensated and compensated density profiles](../figures/compensation_trivial_bulk_example.png)

The compensated amplitude is alpha_0=-28.4253948133-46.6968340029i. The endpoint pump additions are h_L=26.2279471162-47.9655310224i and h_R=-54.6533419295+1.2686970194i; these are additions to the uniform drive epsilon=32. Its uniform bulk winding is -1, with a bulk singular point gap of 0.173666J. The uniform target is exactly stationary and OBC stable. The saved vacuum -> uncompensated relaxation -> smooth compensation ramp protocol reaches it for both T=1000/J and T=10000/J (maximum amplitude errors below 2e-9 after the final hold).

The uncompensated state has frozen winding zero at sites 1–49. Only the last site has a nonzero frozen label; parameters there vary sharply (maximum adjacent coefficient change 1.085J). An isolated boundary label does not establish a bulk topological segment. The smallest-singular-vector centers are near the same right edge (48.15 and 48.52), whereas those of the compensated target lie at opposite ends (48.46 and 2.54). Thus this is evidence of compensation converting an accessible **trivial-bulk** state into a homogeneous topological bulk, rather than a claim of a rigorously classified globally trivial inhomogeneous operator.

We have not excluded other uncompensated stationary branches and cannot assert that this parameter point admits only trivial solutions without compensation. Also, epsilon=38.51 is not an example of trivial uncompensated bulk: its uncompensated attractor already has an extended nonzero frozen-winding segment. The epsilon=32 example makes the narrower, directly tested preparation comparison explicit.

Source: `data/generated/compensation_n50/preparation_30_35/epsilon_32.000000.json`. Reproduce the comparison and 16384-point winding checks with `python src/compensation_trivial_bulk_example.py`; the underlying evolution protocol and tolerances are documented in the preparation scan above.
'''
    path=Path('docs/uniform_topological_segment_stability.md');old=path.read_text(encoding='utf-8')
    old=old.split('## Example: compensation creates an accessible topological bulk at epsilon=32')[0].rstrip()
    path.write_text(old+'\n\n'+text,encoding='utf-8')


if __name__=='__main__':main()
