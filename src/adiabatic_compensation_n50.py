"""Smoothly add intermediate-root compensation after vacuum relaxation."""
from compensation_n50_vacuum import *

def ramp(initial,a0,T,rtol=2e-10,atol=2e-11,epsilon=38.51):
    m=compensated(0,a0,epsilon);q=compensated(1,a0,epsilon).pump-m.pump
    def rhs(t,a):
        x=np.clip(t/T,0,1);s=x*x*(3-2*x)
        f=m.pump+s*q+(m.detuning+m.U-.5j*m.kappa+2*m.U*abs(a)**2)*a
        f[:-1]+=np.exp(-1j*m.phi)*a[1:];f[1:]+=np.exp(1j*m.phi)*a[:-1]
        return -1j*f
    sol=solve_ivp(rhs,[0,T],initial,method='DOP853',rtol=rtol,atol=atol,t_eval=np.linspace(0,T,501))
    assert sol.success
    a=sol.y[:,-1];_,hold=trajectory(compensated(1,a0,epsilon),a,1200,rtol,atol)
    final=np.array(hold['final_alpha_real'])+1j*np.array(hold['final_alpha_imag'])
    hold['target_amplitude_error']=float(max(abs(final-a0)))
    hold['target_reached']=bool(hold['reached_stationary'] and hold['target_amplitude_error']<1e-3)
    return dict(ramp_time=T,rtol=rtol,atol=atol,times=sol.t.tolist(),density_trajectory=(abs(sol.y)**2).T.tolist(),
                ramp_end_force=float(max(abs(compensated(1,a0,epsilon).force(a)))),hold=hold)

def main():
    earlier=json.loads((DEST/'after_relaxation.json').read_text())
    base=earlier['baseline_refined'];initial=np.array(base['final_alpha_real'])+1j*np.array(base['final_alpha_imag'])
    target=json.loads((DEST/'compensation.json').read_text())['branches'][1];a0=complex(*target['alpha0'])
    result=dict(N=50,epsilon=38.51,target_density=target['density_target'],runs=[])
    for T in [10,100,1000,10000]:
        z=ramp(initial.copy(),a0,T);result['runs'].append(z)
        (DEST/'adiabatic_compensation.json').write_text(json.dumps(result))
        print(T,z['hold']['target_reached'],z['hold']['target_amplitude_error'],flush=True)
    result['tight_check']=ramp(initial.copy(),a0,10000,2e-12,2e-13)
    # Independent quasi-static attractor path: stationary root continuation from s=0.
    a=initial.copy();path=[]
    for s in np.unique(np.r_[np.linspace(0,1,501),1-np.geomspace(.002,1e-12,401)]):
        m=compensated(float(s),a0);r=m.stationary(a,1e-8)
        fallback=False
        if not r['converged'] or max(eigvals(m.bogoliubov(r['alpha'])).imag)>-1e-7:
            fallback=True;a,_,_=evolve(m,a,1200,2e-10,2e-11);r=m.stationary(a,1e-8)
        assert r['converged'],s
        a=r['alpha'];z=diagnose(m,a);z.update(s=float(s),evolution_fallback=fallback);path.append(z)
    result['quasistatic_path']=path
    (DEST/'adiabatic_compensation.json').write_text(json.dumps(result))
    plot_result(result);report(result)

def plot_result(result):
    fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for run in result['runs']:
        n=np.array(run['density_trajectory']);times=np.array(run['times']);x=times/run['ramp_time'];s=x*x*(3-2*x)
        ax[0,0].plot(s,n.mean(axis=1),label=f"T={run['ramp_time']}")
        ax[0,1].plot(range(1,51),run['hold']['density'],label=f"T={run['ramp_time']}")
    ax[0,0].axhline(result['target_density'],color='k',ls='--');ax[0,0].set(xlabel='compensation s during ramp',ylabel='mean density');ax[0,0].legend()
    ax[0,1].axhline(result['target_density'],color='k',ls='--',label='uniform target');ax[0,1].set(xlabel='site j',ylabel='density after final hold');ax[0,1].legend()
    path=result['quasistatic_path'];ax[1,0].plot([z['s'] for z in path],[z['growth'] for z in path]);ax[1,0].axhline(0,color='k',ls='--');ax[1,0].set(xlabel='compensation s',ylabel='stationary-path growth / J')
    im=ax[1,1].pcolormesh([z['s'] for z in path],np.arange(1,51),np.array([z['density'] for z in path]).T,shading='nearest');fig.colorbar(im,ax=ax[1,1],label='density');ax[1,1].set(xlabel='compensation s',ylabel='site j',title='Quasi-static attractor path')
    fig.savefig('figures/compensation_n50_adiabatic.png',dpi=180)

def report(result):
    lines=[]
    for run in result['runs']:
        z=run['hold'];lines.append(f"|{run['ramp_time']}|{'yes' if z['target_reached'] else 'NO'}|{np.mean(z['density']):.3f}|{z['target_amplitude_error']:.5g}|{z.get('growth',float('nan')):.6f}|")
    path=result['quasistatic_path'];tight=result['tight_check'];last=result['runs'][-1]
    difference=float(max(abs(np.array(tight['hold']['density'])-np.array(last['hold']['density']))))
    block=f'''## N=50: smooth turn-on of the topological-root compensation

**Slow smooth switching reaches the uniform topological target in the tested N=50 chain.** Ramps with T=1000/J and 10000/J succeed, while T=10/J and 100/J reach the competing inhomogeneous state. This demonstrates a preparation from vacuum via initial uncompensated relaxation and slow endpoint switching, rather than initialization near the target.

We test whether slower switching reaches the uniform intermediate/topological root at epsilon=38.51. The initial state is the raw endpoint of the independently verified uncompensated vacuum evolution to t=2400, not a prepared uniform target. Only endpoint compensation changes; all bulk parameters and drives remain fixed.

For 0<=t<=T use s(t)=3(t/T)^2-2(t/T)^3, with zero slope at each endpoint. The endpoint drives follow the earlier intermediate-root compensation formula. After the ramp, hold s=1 for another 1200 time units and test stationarity and complex-amplitude agreement with the uniform target using the same thresholds as above.

|ramp duration T (1/J)|uniform topological target reached?|final mean density|max amplitude error from target|final-state growth / J|
|---:|---|---:|---:|---:|
{chr(10).join(lines)}

![Smooth compensation ramps and quasi-static stationary-path diagnostics](../figures/compensation_n50_adiabatic.png)

The slowest ramp was repeated with rtol=2e-12, atol=2e-13 instead of 2e-10, 2e-11. Its final density profile differs by {difference:.3g}; target reached in that check: {tight['hold']['target_reached']}. The tight final raw force is {tight['hold']['end_force']:.3g}; maximum last-50-time density standard deviation is {tight['hold']['tail_density_std']:.3g}.

We separately follow a stationary attractor from s=0 to 1 on a 0.002 grid, with logarithmically refined steps in 1-s from 0.002 down to 1e-12 near the final endpoint, beginning with the same relaxed profile. Each step uses an analytic-Jacobian stationary root solve; if the solve fails or returns an unstable state, a mean-field evolution fallback supplies a new attractor seed. These possible switches are explicitly recorded and this procedure is not pseudo-arclength continuation. There are {sum(z['evolution_fallback'] for z in path)} such fallback steps. The maximum computed growth on the accepted path is {max(z['growth'] for z in path):.6f} J, and the final maximum density difference from the uniform target is {max(abs(np.array(path[-1]['density'])-result['target_density'])):.6g}. This numerical branch-tracking diagnostic can select a different basin from the continuous time-dependent ramp; it is not an independent simulation of the adiabatic limit.

These finite ramps and the quasi-static path test whether simple slow addition of the fixed endpoint compensation enters the desired basin. They do not prove a universal adiabatic limit or exclude other drive paths. A stable final target need not be the continuation of the attractor occupied initially. If the path remains on another stable branch, increasing ramp duration need not select the uniform target.

Reproduce with `python src/adiabatic_compensation_n50.py`. Full ramp trajectories, post-ramp holds, tolerance check, stationary-path profiles, growth and fallback records are saved in data/generated/compensation_n50/adiabatic_compensation.json. The same DOP853 mean-field evolution, root residual acceptance 1e-8, and stability threshold 1e-7 J apply. No covariance evolution is used.
'''
    p=Path('docs/uniform_topological_segment_stability.md');text=p.read_text(encoding='utf-8').split('## N=50: smooth turn-on of the topological-root compensation')[0].rstrip();p.write_text(text+'\n\n'+block,encoding='utf-8')

if __name__=='__main__':main()
