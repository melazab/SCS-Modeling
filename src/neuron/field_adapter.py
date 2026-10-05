"""Small, validated per-compartment field files; no FEM mesh in NEURON workers."""
from pathlib import Path
import hashlib
import json
import numpy as np
from stimulation import combine_basis


def trajectory_positions(polyline_mm, s_um):
    p = np.asarray(polyline_mm,float)
    s = np.asarray(s_um,float)
    if p.ndim != 2 or p.shape[1] != 3 or len(p) < 2 or not np.isfinite(p).all():
        raise ValueError('Trajectory must contain at least two finite 3D points in mm')
    lengths = np.linalg.norm(np.diff(p,axis=0),axis=1)
    if np.any(lengths <= 0) or not np.isfinite(s).all() or s.ndim != 1 or np.any(s < 0):
        raise ValueError('Invalid trajectory or compartment arc lengths')
    arc = np.r_[0,np.cumsum(lengths)]
    if s.max(initial=0)*.001 > arc[-1] + 1e-9:
        raise ValueError('Axon exceeds trajectory length; endpoints will not be clamped')
    return np.column_stack([np.interp(s*.001,arc,p[:,i]) for i in range(3)])


def geometry_hash(kind, s_um, xyz_mm):
    h = hashlib.sha256(kind.encode())
    for array in (s_um,xyz_mm):
        a = np.asarray(array,dtype='<f8')
        h.update(str(a.shape).encode())
        h.update(a.tobytes())
    return h.hexdigest()


def load_transfer(path, kind, s_um, xyz_mm, contact_weights):
    with np.load(path,allow_pickle=False) as d:
        if str(d['units_coords']) != 'mm' or str(d['units_phi']) != 'V/A':
            raise ValueError('Unexpected field units')
        if str(d['geometry_hash']) != geometry_hash(kind,s_um,xyz_mm):
            raise ValueError('Sampled field belongs to different axon coordinates')
        if not np.array_equal(d['s_um'],s_um) or not np.array_equal(d['xyz_mm'],xyz_mm):
            raise ValueError('Compartment coordinates differ')
        ids = [int(i) for i in d['contact_ids']]
        if len(set(ids)) != len(ids) or set(contact_weights) - set(ids):
            raise ValueError('Unknown or duplicate contact IDs')
        field = combine_basis(d['phi_V_A'],[contact_weights.get(i,0.) for i in ids])
        meta = json.loads(str(d['metadata_json']))
        return field,meta
