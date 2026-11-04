"""Refine existing stationary-path diagnostic without rerunning time ramps."""
from adiabatic_compensation_n50 import *

def main():
    p=DEST/'adiabatic_compensation.json';r=json.loads(p.read_text())
    target=json.loads((DEST/'compensation.json').read_text())['branches'][1];a0=complex(*target['alpha0'])
    path=[z for z in r['quasistatic_path'] if z['s']<=.998]
    a=np.array(path[-1]['alpha_real'])+1j*np.array(path[-1]['alpha_imag'])
    for s in np.r_[1-np.geomspace(.002,1e-12,401)[1:],1.]:
        m=compensated(float(s),a0);z=m.stationary(a,1e-8);fallback=False
        if not z['converged'] or max(eigvals(m.bogoliubov(z['alpha'])).imag)>-1e-7:
            fallback=True;a,_,_=evolve(m,a,1200,2e-10,2e-11);z=m.stationary(a,1e-8)
        assert z['converged'];a=z['alpha'];rec=diagnose(m,a);rec.update(s=float(s),evolution_fallback=fallback);path.append(rec)
    r['quasistatic_path']=path;p.write_text(json.dumps(r));plot_result(r);report(r)
    print('end target amplitude error',max(abs(a-a0)),flush=True)

if __name__=='__main__':main()
