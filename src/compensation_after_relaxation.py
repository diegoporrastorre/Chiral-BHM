"""Vacuum -> uncompensated steady state -> fixed endpoint compensation."""
from compensation_n50_vacuum import *

def compare(z,a0):
    a=np.array(z['final_alpha_real'])+1j*np.array(z['final_alpha_imag'])
    z['target_amplitude_error']=float(max(abs(a-a0)))
    z['target_density_error']=float(max(abs(abs(a)**2-abs(a0)**2)))
    z['target_reached']=bool(z['reached_stationary'] and z['target_amplitude_error']<1e-3)

def main():
    previous=json.loads((DEST/'compensation.json').read_text())
    m=Chain(N=50,drive=38.51)
    a,base=trajectory(m,np.zeros(50,complex))
    aa,base_refined=trajectory(m,np.zeros(50,complex),2400,2e-12,2e-13)
    assert base['reached_stationary'] and base_refined['reached_stationary']
    result=dict(N=50,epsilon=38.51,baseline=base,baseline_refined=base_refined,branches=[])
    for b in previous['branches']:
        a0=complex(*b['alpha0']);m=compensated(1,a0)
        # Copy raw integration endpoints, not Newton-polished substitutes.
        _,z=trajectory(m,a.copy())
        _,zz=trajectory(m,aa.copy(),2400,2e-12,2e-13)
        compare(z,a0);compare(zz,a0)
        rec=dict(label=b['label'],target_density=b['density_target'],endpoint_pumps=b['endpoint_pumps'],
                 after_switch=z,after_switch_refined=zz)
        result['branches'].append(rec)
        (DEST/'after_relaxation.json').write_text(json.dumps(result))
        print(b['label'],z['target_reached'],zz['target_reached'],z['target_amplitude_error'],flush=True)
    plot_result(result);report(result)

def plot_result(result):
    fig,ax=plt.subplots(2,3,figsize=(14,8),layout='constrained')
    for col,b in enumerate(result['branches']):
        z=b['after_switch'];t=np.array(z['times']);n=np.array(z['density_trajectory'])
        for j in [0,9,24,39,49]:ax[0,col].plot(t,n[:,j],label=f'j={j+1}')
        ax[0,col].axhline(b['target_density'],color='k',ls='--',label='uniform target')
        ax[0,col].set(xlabel='time after compensation switch Jt',ylabel='density',title=b['label']+' compensation',xlim=(0,150));ax[0,col].legend(fontsize=7,ncol=2)
        ax[1,col].plot(range(1,51),result['baseline']['density'],color='.6',ls=':',label='before compensation')
        ax[1,col].plot(range(1,51),z['density'],label='after switch, t=1200')
        ax[1,col].plot(range(1,51),b['after_switch_refined']['density'],ls='--',label='t=2400, tighter tolerances')
        ax[1,col].axhline(b['target_density'],color='k',ls='--',label='uniform target')
        ax[1,col].set(xlabel='site j',ylabel='density');ax[1,col].legend(fontsize=7)
    fig.savefig('figures/compensation_n50_after_relaxation.png',dpi=180)

def report(result):
    lines=[]
    for b in result['branches']:
        z=b['after_switch'];zz=b['after_switch_refined']
        lines.append(f"|{b['label']}|{b['target_density']:.3f}|{'yes' if z['target_reached'] else 'NO'}|{'yes' if zz['target_reached'] else 'NO'}|{z['target_amplitude_error']:.4g}|{z.get('growth',float('nan')):.6f}|")
    mid=result['branches'][1];z=mid['after_switch'];zz=mid['after_switch_refined']
    original=json.loads((DEST/'compensation.json').read_text())['branches'][1]['vacuum']
    difference=float(max(abs(np.array(z['density'])-np.array(original['density']))))
    text=f'''## N=50: add compensation after uncompensated relaxation

At epsilon=38.51 we now test the proposed two-stage preparation, rather than applying compensation directly to vacuum. All model parameters remain as in the preceding N=50 calculation.

1. Start at alpha_j=0 and evolve the uniformly driven **uncompensated** OBC chain to t=1200.
2. Verify that this raw dynamical endpoint is stationary, then use precisely that same complex profile as the initial condition in three separate experiments.
3. Abruptly add the fixed low-, intermediate-, or high-root endpoint compensation, keep all drives constant thereafter, and evolve for another 1200 time units.
4. Repeat the entire two-stage experiment independently with 2400 time units per stage and tighter tolerances.

No target state is seeded and no root-polished state replaces the dynamical endpoint at the switch. Root solving is only a final diagnostic. Compensation changes only the endpoint drive amplitudes/phases; it does not change detuning, hopping, loss, or interior drive.

The first stage reaches a stable inhomogeneous stationary state with mean density {np.mean(result['baseline']['density']):.3f}, growth {result['baseline']['growth']:.6f} J and raw endpoint force {result['baseline']['end_force']:.3g}. The two preparation runs agree in their final site densities to {max(abs(np.array(result['baseline']['density'])-np.array(result['baseline_refined']['density']))):.3g}.

|compensation target|target uniform density|uniform target reached, 1200 per stage|target reached, tighter 2400 per stage|max complex amplitude error from target, 1200|final-state growth / J|
|---|---:|---|---|---:|---:|
{chr(10).join(lines)}

![Switching on endpoint compensation after reaching the uncompensated stationary profile](../figures/compensation_n50_after_relaxation.png)

The intermediate compensation ends with mean density {np.mean(z['density']):.3f}, range {min(z['density']):.3f}-{max(z['density']):.3f}. Its density difference from the inhomogeneous attractor obtained by applying intermediate compensation directly to vacuum in the preceding experiment is {difference:.3g}. The tighter second-stage run has raw final force {zz['end_force']:.3g}, maximum amplitude distance to the root-polished endpoint {zz['root_distance']:.3g}, and maximum density standard deviation over its last 50 time units {zz['tail_density_std']:.3g}.

The results distinguish the basin reached by this protocol from local linear stability of the uniform targets. Failure of the intermediate target to be reached here would not contradict its verified local stability, nor rule out other preparation paths. Successful targets are identified by complex amplitude, not density alone: max_j|alpha_j-alpha0|<1e-3, in addition to the stationarity checks described above.

The base runs use DOP853 with rtol=2e-10, atol=2e-11; the independent longer runs use rtol=2e-12, atol=2e-13. Raw trajectories, initial preparation, final-state diagnostics and target errors are stored in data/generated/compensation_n50/after_relaxation.json. Reproduce with `python src/compensation_after_relaxation.py` from the project root. This session used C:/Users/diego/miniconda3/envs/topo/python.exe. No Gaussian covariance equations are evolved.
'''
    p=Path('docs/uniform_topological_segment_stability.md');base=p.read_text(encoding='utf-8').split('## N=50: add compensation after uncompensated relaxation')[0].rstrip();p.write_text(base+'\n\n'+text,encoding='utf-8')

if __name__=='__main__':main()
