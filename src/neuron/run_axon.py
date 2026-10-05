"""Run one representative axon against an analytic field or cached Elmer field."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys
import numpy as np
from axons import Axon
from field_adapter import trajectory_positions,load_transfer
from simulation import trial,threshold,crossings
from stimulation import Pulse,point_source_basis,combine_basis

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def save_trace(path,result):
    np.savez_compressed(path,t_ms=result['t_ms'],v_mV=result['v_mV'],s_um=result['s_um'])


def run(args):
    out=Path(args.out).resolve();out.mkdir(parents=True,exist_ok=True)
    axon=Axon(args.kind,args.length_um,args.refinement,args.temperature)
    pulse=Pulse(width_ms=args.width_ms)
    try:
        # Intracellular propagation is an independent capability check, with
        # the original small-fibre 1-ms pulse rather than the threshold pulse.
        intra=trial(axon,intracellular_nA=2. if args.kind=='abeta' else .2,
                    pulse=Pulse(width_ms=1),dt_ms=args.dt_ms)
        save_trace(out/'propagation.npz',intra)
        t1=crossings(intra['t_ms'],intra['v_mV'][1],5)
        t2=crossings(intra['t_ms'],intra['v_mV'][7],5)
        if not t1.size or not t2.size or t2[0] <= t1[0] or intra['spontaneous']:
            raise RuntimeError('Intracellular propagation validation failed')
        velocity=float((intra['s_um'][7]-intra['s_um'][1])/(t2[0]-t1[0])*.001)
        if args.solution:
            if not args.trajectory:
                raise ValueError('--trajectory is required with --solution')
            trajectory=json.loads(Path(args.trajectory).read_text())
            if trajectory.get('units') != 'mm':
                raise ValueError('Trajectory must explicitly specify mm')
            if sum(np.linalg.norm(np.diff(trajectory['points'],axis=0),axis=1))*1000 < axon.length_um:
                raise ValueError('Trajectory must cover the full constructed axon length')
            xyz=trajectory_positions(trajectory['points'],axon.s_um)
            request=out/'sample_request.npz'
            np.savez_compressed(request,kind=args.kind,s_um=axon.s_um,xyz_mm=xyz,
                                units_coords='mm',trajectory_json=json.dumps(trajectory))
            subprocess.run([str(ROOT/'fem/.venv/bin/python'),str(HERE/'sample_fem.py'),
                            str(Path(args.solution).resolve()),str(request),str(out/'sampled_field.npz')],check=True)
            weights={int(k):float(v) for k,v in json.loads(args.weights).items()}
            transfer,metadata=load_transfer(out/'sampled_field.npz',args.kind,axon.s_um,xyz,weights)
            metadata['contact_weights']=weights
        else:
            xyz=np.column_stack((axon.s_um,np.zeros((len(axon.s_um),2))))
            electrodes=np.array([[.4*axon.length_um,args.distance_um,0],[.6*axon.length_um,args.distance_um,0]])
            transfer=combine_basis(point_source_basis(xyz,electrodes),[-1,1])
            metadata=dict(source='analytic bipolar point sources',sigma_S_m=.2,
                          electrodes_um=electrodes.tolist(),weights=[-1,1])
        summary,result=threshold(axon,transfer,pulse,args.dt_ms,args.max_ma)
        save_trace(out/'threshold_trace.npz',result)
        np.savez_compressed(out/'transfer.npz',s_um=axon.s_um,transfer_V_A=transfer)
        baseline=trial(axon,transfer,0,pulse,args.dt_ms)
        uniform=trial(axon,np.full(len(axon.s_um),1000.),1.,pulse,args.dt_ms)
        gauge_error=float(np.max(np.abs(uniform['v_mV']-baseline['v_mV'])))
        if baseline['activated'] or uniform['activated'] or gauge_error > 1e-4:
            raise RuntimeError(f'Zero/uniform field control failed: {gauge_error:g} mV')
        import neuron
        summary.update(kind=args.kind,diameter_um=axon.diameter_um,length_um=axon.length_um,
            temperature_C=axon.temperature_C,compartments=len(axon.s_um),dt_ms=args.dt_ms,
            refinement=args.refinement,pulse=dict(start_ms=pulse.start_ms,width_ms=pulse.width_ms,
                gap_ms=pulse.gap_ms,phases=[1,-1],units='ms',description='symmetric charge-balanced biphasic'),
            conduction_velocity_m_s=velocity,uniform_field_error_mV=gauge_error,
            arrivals=result['arrivals'],field=metadata,neuron_version=neuron.__version__,
            sources=json.loads((HERE/'vendor/sources.json').read_text()),
            code_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in HERE.glob('*.py')})
        (out/'report.json').write_text(json.dumps(summary,indent=2)+'\n')
        print(json.dumps({k:summary[k] for k in ('kind','status','threshold_mA','conduction_velocity_m_s','uniform_field_error_mV')},indent=2),flush=True)
    finally:
        axon.close()


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--kind',choices=['abeta','adelta','c'],required=True)
    p.add_argument('--out',required=True)
    p.add_argument('--length-um',type=float,default=10000)
    p.add_argument('--refinement',type=int,choices=[1,2],default=1)
    p.add_argument('--temperature',type=float,default=37.)
    p.add_argument('--width-ms',type=float,default=.3)
    p.add_argument('--dt-ms',type=float,default=.005)
    p.add_argument('--max-ma',type=float,default=100.)
    p.add_argument('--distance-um',type=float,default=500.)
    p.add_argument('--solution');p.add_argument('--trajectory')
    p.add_argument('--weights',default='{"1":-1,"3":1}',help='Contact ID to relative current; amplitude multiplies these weights')
    return p


if __name__ == '__main__':
    run(parser().parse_args())
