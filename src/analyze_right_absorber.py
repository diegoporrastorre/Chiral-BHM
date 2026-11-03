"""Precision-refined modal reflection and plots for the right-layer search."""
import json
from pathlib import Path
import numpy as np
from bdg_right_reflection import reflection_decimal

ROOT=Path(__file__).resolve().parents[1]


def run():
    omega=[0,.001,.003,.01,.03,.05,.1,.2,.3,.5,.75,1,1.5,2,3]
    cases=[(0,'quintic',0),(80,'quintic',0),(90,'quintic',0),(96,'quintic',0),(80,'cap4',10),(90,'cap2',10)]
    output=dict(kappa=2.1,frequency_units='J',normalization='unit Euclidean norm bulk mode vectors; not flux',rows=[])
    for l,f,e in cases:
        coarse=reflection_decimal(l,f,e,2.1,omega,100)
        fine=reflection_decimal(l,f,e,2.1,omega,160)
        error=max(abs(a['abs_r']-b['abs_r'])/max(1.,b['abs_r']) for a,b in zip(coarse,fine))
        assert error<1e-12
        output['rows'].append(dict(length=l,family=f,eta=e,precision_comparison=error,response=fine))
        print('REFLECTION',l,f,e,[(r['omega'],r['abs_r']) for r in fine if r['omega'] in [0,.5,1]],flush=True)
    transparent=reflection_decimal(4,'loss',0,2.1,[0,.5,1,2],100,True)
    assert max(r['abs_r'] for r in transparent)<1e-85
    output['transparent_reference']=transparent
    # Independent boundary-mode quantization check at a known sharp N=8 root.
    check=reflection_decimal(0,'quintic',0,3.300416850278589,[0],100,n=8)[0]
    mismatch=abs(check['abs_r']/check['unit_loop_reflection_scale']-1)
    assert mismatch<1e-10
    output['sharp8_loop_check_relative_error']=mismatch
    (ROOT/'data/generated/bdg_right_reflection.json').write_text(json.dumps(output,indent=2))
    return output


def plot(response):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    search=json.loads((ROOT/'data/generated/bdg_right_absorber.json').read_text())
    fig,axs=plt.subplots(2,2,figsize=(11,7.5),layout='constrained')
    ax=axs[0,0]
    for family in ['linear','cubic','quintic','tanh2','tanh4']:
        rows=sorted([r for r in search['rows'] if r['family']==family and r['eta']==0 and r['length']>0],key=lambda r:r['length'])
        ax.plot([r['length'] for r in rows],[r['critical_kappa'] for r in rows],'o-',label=family)
    ax.axhline(2,color='k',ls='--',label='Original bulk threshold')
    ax.set(xlabel='Right-layer length L (sites)',ylabel=r'$\kappa_c/J$',title='Sharp left edge; exact uniform core of 100-L sites')
    ax.legend(fontsize=7,ncol=2)
    ax=axs[0,1]
    for l,f,e in [(80,'quintic',0),(90,'quintic',0),(90,'cap2',10)]:
        from bdg_right_absorber import profiles
        s,loss=profiles(100,l,f,e)
        if f=='quintic': ax.plot(np.arange(1,101),[float(v*v) for v in s],label=f'Pairing/detuning, L={l}')
        else: ax.plot(np.arange(1,101),[float(v)/e for v in loss],ls='--',label='Extra amplitude loss / 10, L=90')
    ax.set(xlabel='Site j',ylabel='Profile factor',title='Only the right side is modified')
    ax.legend(fontsize=8)
    ax=axs[1,0]
    for row in response['rows']:
        if row['length']==96: continue
        rr=row['response']
        label='Sharp boundary' if row['length']==0 else f"L={row['length']} {row['family']}"
        ax.semilogy([r['omega'] for r in rr],[r['abs_r'] for r in rr],'o-',ms=3,label=label)
    ax.set(xlabel=r'$\omega/J$',ylabel='Returning / incident modal amplitude',title=r'Sampled right reflection at $\kappa=2.1J$')
    ax.legend(fontsize=7)
    ax=axs[1,1]
    for l in [80,90,96]:
        row=next(r for r in response['rows'] if r['length']==l and r['family']=='quintic')
        rr=row['response']
        ax.semilogy([r['omega'] for r in rr],[r['abs_r']/r['unit_loop_reflection_scale'] for r in rr],'o-',ms=3,label=f'L={l}, uniform core {100-l}')
    ax.axhline(1,color='k',ls='--')
    ax.set(xlabel=r'$\omega/J$',ylabel='Round-trip amplitude magnitude',title='Feedback magnitude (phase also matters)')
    ax.legend(fontsize=8)
    fig.suptitle(r'$N=100$, bulk $\Delta=0.1J$, $g_s=g_c=J=1$, $\phi=\pi/2$: right-only absorbers')
    fig.savefig(ROOT/'figures/bdg_right_absorber.png',dpi=180)
    plt.close(fig)


if __name__=='__main__': plot(run())
