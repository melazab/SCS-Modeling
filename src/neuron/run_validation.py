"""Reproduce three-fibre analytic validation in isolated NEURON processes."""
from pathlib import Path
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
import subprocess
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
KINDS=('abeta','adelta','c')
CASES=(('base',()),('time',('--dt-ms','.0025')),('space',('--refinement','2')))


def run_case(out,kind,name,options):
    directory=out/kind/name
    completed=subprocess.run(
        [sys.executable,str(HERE/'run_axon.py'),'--kind',kind,
         '--out',str(directory),*options],
        check=False,capture_output=True,text=True)
    if completed.returncode:
        raise RuntimeError(
            f'{kind}/{name} failed (exit {completed.returncode})\n'
            f'{completed.stdout}{completed.stderr}')
    return kind,name,json.loads((directory/'report.json').read_text()),completed.stdout


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',default=str(HERE/'out/validation'))
    default_jobs=min(len(KINDS)*len(CASES),max(1,(os.cpu_count() or 1)//2))
    p.add_argument('--jobs',type=int,default=default_jobs,
                   help=f'independent NEURON processes (default: {default_jobs})')
    a=p.parse_args();out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=True)
    if a.jobs < 1:
        p.error('--jobs must be at least 1')
    reports={kind:{} for kind in KINDS}
    jobs=[(kind,name,options) for kind in KINDS for name,options in CASES]
    print(f'Running {len(jobs)} validation cases with {min(a.jobs,len(jobs))} processes',flush=True)
    with ThreadPoolExecutor(max_workers=min(a.jobs,len(jobs))) as executor:
        futures={executor.submit(run_case,out,*job):job for job in jobs}
        for future in as_completed(futures):
            kind,name,report,stdout=future.result()
            reports[kind][name]=report
            print(f'completed {kind}/{name}',flush=True)
            if stdout.strip():
                print(stdout.strip(),flush=True)
    summary={}
    for kind in KINDS:
        kind_reports=reports[kind]
        base=kind_reports['base']['threshold_mA']
        if base is None:raise RuntimeError(f'{kind}: no analytic threshold found')
        errors={name:abs(kind_reports[name]['threshold_mA']/base-1)
                for name,_ in CASES if name!='base'}
        # Engineering convergence gate, not physiological validation.
        if any(e>.05 for e in errors.values()):
            raise RuntimeError(f'{kind}: >5% threshold change on refinement: {errors}')
        summary[kind]=dict(threshold_mA=base,velocity_m_s=kind_reports['base']['conduction_velocity_m_s'],
                           relative_refinement_changes=errors)
    plot(out)
    temporary=out/'summary.json.tmp'
    temporary.write_text(json.dumps(summary,indent=2)+'\n')
    temporary.replace(out/'summary.json')


def plot(out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(3,2,figsize=(10,8),layout='constrained')
    for row,kind in enumerate(KINDS):
        directory=out/kind/'base'
        for col,name in enumerate(('propagation','threshold_trace')):
            with np.load(directory/f'{name}.npz') as d:
                for i in (1,4,7):axes[row,col].plot(d['t_ms'],d['v_mV'][i],label=f"{d['s_um'][i]/1000:.1f} mm")
            axes[row,col].set(xlabel='Time (ms)',ylabel='Membrane potential (mV)',
                title=f'{kind.upper()} — '+('intracellular propagation' if col==0 else 'analytic field threshold'))
            axes[row,col].legend(fontsize=8)
    fig.suptitle('Three representative axons — engineering validation, not SCS recruitment predictions')
    fig.savefig(out/'validation.png',dpi=180)
    fig.savefig(out/'validation.pdf')
    plt.close(fig)


if __name__ == '__main__':main()
