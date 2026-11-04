"""Append or replace the fine-scan section in the existing stability report."""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/generated/homogeneous_obc_n100'

def main():
    s=json.loads((DATA/'fine_scan_summary.json').read_text())
    rows=json.loads((DATA/'fine_scan_35_45.json').read_text())
    up=rows[:1001]
    def longest(r):
        best=run=0
        for nu in r.get('nu',[]):
            run=run+1 if nu else 0;best=max(best,run)
        return best
    extended=[r for r in up if longest(r)>=5]
    table=[]
    for eps in [38.45,38.46,38.47,38.48,38.49,38.50,38.51,38.52,38.53,38.54,38.55]:
        r=next(r for r in up if r['epsilon']==eps)
        n=r['density']
        table.append(f"|{eps:.2f}|{sum(n)/100:.2f}|{longest(r)}|{r['growth']:.6f}|{r['smin']:.6g}|")
    w=s['worst_uniform'];different=s['different_drives']
    section=f'''## Finer scan: epsilon spacing 0.01 across 35-45

Added 2026-11-04. Parameters, homogeneous drive/detuning, N=100, OBC, V=0, and all sign and topology conventions remain as above. The full interval is now scanned at **0.01 spacing**, rather than 0.25. Each direction contains 1001 points, including endpoints 35 and 45 as continuation anchors; the strict interior contains 999 points per direction. Both directions were recalculated independently from their respective coarse-scan endpoint solutions.

### Convergence and full-chain stability

Across {s['points']} stationary attempts there are {s['failed']} failed solves; maximum stationary residual is {s['max_residual']:.3g}. The largest full-chain growth rate is {s['max_growth']:.6f} J. Stability classification uses tolerance 1e-7 J. Root tolerance is 1e-10, acceptance max|f|<1e-8; failed or unstable roots trigger mean-field evolution to t=1200 with rtol=2e-10, atol=2e-11 before another solve. No covariance equations are evolved.

The largest up/down density difference is {s['max_direction_density_difference']:.6g}. There are {len(different)} sampled drives with maximum site-density difference above 0.01 between directions. {'Their drives are '+', '.join(f'{x:.2f}' for x in different)+'.' if len(different)<=30 and different else ''} These are finite-step continuation comparisons, not an exhaustive branch enumeration.

![Fine density and local winding maps in both directions](../figures/fine_scan_35_45_maps.png)

### Resolving the narrow topological crossover

The up-sweep's most extended slowly varying candidate occurs at epsilon={s['best_slow_domain_epsilon']:.2f}: {len(s['best_slow_sites'])} sites have both nonzero frozen winding and adjacent-bond variation below 0.05J. This uses the same descriptive variation criterion as the earlier report, not a controlled error bound near local gap closures. The first and last drives with a contiguous winding segment of at least five sites are {extended[0]['epsilon']:.2f} and {extended[-1]['epsilon']:.2f}; this span alone does not imply a single homogeneous phase or certify slow variation throughout it.

|epsilon|mean density|longest nonzero-winding segment (sites)|full-chain growth / J|smallest singular value / J|
|---:|---:|---:|---:|---:|
{chr(10).join(table)}

![Density profiles at closely spaced drives around the largest slowly varying segment; amber dots denote local nonzero winding](../figures/fine_scan_35_45_profiles.png)

The profiles are full nonlinear stationary backgrounds; the uniform chains below are fixed-coefficient quadratic comparisons. The dots include short boundary signals and therefore retain the earlier local-topology caveats.

### Does the uniform-OBC stability conclusion survive refinement?

At each strict-interior up-sweep point we select the longest contiguous candidate segment and repeat its minimum, median, and maximum density coefficients uniformly. Both its own length and length 100 are tested. The {s['uniform_cases']} resulting spectra contain {s['uniform_unstable']} growing cases above 1e-7 J. The least negative rate is {w['growth']:.6f} J, at epsilon={w['epsilon']:.2f}, density {w['density']:.3f}, length {w['N']}, using the {w['sample']} density sample. This is a representative-coefficient search, not proof for all densities between those samples. Only the longest segment is retested on the fine grid; the earlier coarser calculation included all segments.

![Uniformized segment stability on the fine grid, and full inhomogeneous-chain stability](../figures/fine_scan_35_45_stability.png)

For the least stable uniform example, independent Decimal characteristic-polynomial/Routh calculations at 100 and 160 digits both find {s['high_precision'][0]['counts']['0']} growing generator eigenvalues to the right of zero. Counts also agree at shifts -0.1 and +0.1. This checks its stable/unstable label despite non-Hermitian eigenvalue conditioning; it does not certify every spectrum with interval arithmetic. Previously established high-precision checks for epsilon=38.50 still apply.

The fine calculation therefore tests whether the earlier stable uniform-OBC examples persist under closer sampling. It does not change the distinction between a uniform quadratic operator and a newly solved nonlinear chain after cutting boundaries, nor does it imply that stability suppresses transient or noise amplification.

### Reproduce the additional calculation

```powershell
python src/fine_scan_topological_segments.py
python src/report_fine_scan.py
```

The scan is resumable from saved checkpoints. Session interpreter: C:/Users/diego/miniconda3/envs/topo/python.exe. Saved records are fine_scan_35_45.json (complex amplitudes, densities, convergence, full-chain growth, winding, SVD and variation), fine_uniform_segment_tests.json (sampled uniform coefficients, growth, residuals and conditioning), and fine_scan_summary.json (comparison and precision checks), all under data/generated/homogeneous_obc_n100/. The coarse figures and results above are retained for comparison; use this section for the refined scan.
'''
    p=ROOT/'docs/uniform_topological_segment_stability.md'
    text=p.read_text(encoding='utf-8').split('## Finer scan: epsilon spacing 0.01 across 35-45')[0].rstrip()
    p.write_text(text+'\n\n'+section,encoding='utf-8')

if __name__=='__main__':main()
