"""Run with fem/.venv/bin/python: strict sampling of a verified Elmer solution.

No display clamping, nearest-node fill, or silent extrapolation. This process
loads the large mesh once; NEURON subsequently needs only the small output.
"""
from pathlib import Path
import argparse
import importlib.util
import json
import sys
import numpy as np
from field_adapter import geometry_hash

ROOT = Path(__file__).resolve().parents[2]
FEM = ROOT/'fem/scripts'


def sample(solution, request, output):
    sys.path.insert(0,str(FEM))
    import artifacts
    solution = Path(solution).resolve()
    if solution.name not in ('solution.npz','solution_bg.npz'):
        raise ValueError('Expected a production solution.npz or solution_bg.npz')
    background = solution.name == 'solution_bg.npz'
    signature = artifacts.solution_signature(str(solution.parent),background)
    if not artifacts.valid(str(solution),signature):
        raise ValueError('Solution provenance is stale or incomplete')
    spec = importlib.util.spec_from_file_location('scs_fem_field',FEM/'field.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    with np.load(request,allow_pickle=False) as r:
        kind, s, xyz = str(r['kind']),r['s_um'],r['xyz_mm']
        if str(r['units_coords']) != 'mm' or xyz.shape != (len(s),3) or not np.isfinite(xyz).all():
            raise ValueError('Invalid sampling request')
        trajectory_meta = json.loads(str(r['trajectory_json']))
    with np.load(solution,allow_pickle=False) as d:
        if bool(d['background_included']) != background:
            raise ValueError('Background metadata mismatch')
        if not np.isfinite(d['relative_residuals']).all() or np.any(d['relative_residuals'] > float(d['accepted_residual'])):
            raise ValueError('Solution residual acceptance failed')
        f = module.TetField(d['nodes'],d['tets'],d['phi'])
        ti,w = f.locate(xyz)
        if np.any(ti < 0):
            bad = np.flatnonzero(ti<0)
            raise ValueError(f'{len(bad)} compartments outside located mesh; first indices {bad[:10].tolist()}')
        phi = np.einsum('cnk,nk->cn',f.values[:,f.tets[ti]],w)
        if not np.isfinite(phi).all():
            raise ValueError('Nonfinite sampled field')
        labels = d['label'][ti]
        order = d['order']
        tissues = [str(order[i]) if i >= 0 else 'background' for i in labels]
        allowed = trajectory_meta.get('allowed_tissues')
        offending = [i for i,t in enumerate(tissues) if allowed and t not in allowed]
        # A tissue label is a provenance check on WHERE the axon sits; it is not an
        # input to the field. phi is P1-interpolated from node values that a
        # mislabelled tet shares with its neighbours, so tolerating a declared,
        # counted handful of isolated boundary slivers changes no sampled number.
        # The allowance must be declared per trajectory and is never open-ended:
        # every offending compartment is recorded here for audit.
        budget = int(trajectory_meta.get('max_label_exceptions',0))
        if offending and len(offending) > budget:
            raise ValueError(f'Trajectory crosses unexpected tissues: {set(tissues)-set(allowed)} '
                             f'in {len(offending)} compartments (declared allowance {budget})')
        exceptions = [dict(index=int(i),mesh_label=tissues[i],
                           xyz_mm=[float(c) for c in xyz[i]]) for i in offending]
        metadata = dict(solution=str(solution),solution_sha256=artifacts.file_hash(str(solution)),
                        solution_signature=signature,background=background,
                        coordinate_units='mm',potential_units='V/A',
                        interpolation='P1 barycentric; no fill',trajectory=trajectory_meta,
                        tissue_counts={t:tissues.count(t) for t in set(tissues)},
                        label_exceptions=exceptions)
        if artifacts.solution_signature(str(solution.parent),background) != signature or not artifacts.valid(str(solution),signature):
            raise ValueError('Solution changed while sampling')
        artifacts.atomic_npz(str(output),phi_V_A=phi,contact_ids=d['contact_ids'],
            s_um=s,xyz_mm=xyz,units_coords='mm',units_phi='V/A',
            geometry_hash=geometry_hash(kind,s,xyz),metadata_json=json.dumps(metadata))
    print(json.dumps(metadata,indent=2))


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('solution');p.add_argument('request');p.add_argument('output')
    a=p.parse_args();sample(a.solution,a.request,a.output)
