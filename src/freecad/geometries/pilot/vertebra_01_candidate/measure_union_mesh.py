"""Diagnostic mesh union; this does not repair a failed native CAD Boolean."""
import json,time
from pathlib import Path
import numpy as np,trimesh
HERE=Path(__file__).resolve().parent
started=time.perf_counter()
a=trimesh.load_mesh(HERE/'body_candidate.stl'); d=json.loads((HERE/'posterior_working.json').read_text()); b=trimesh.Trimesh(d['vertices'],d['faces'],process=False)
c=trimesh.boolean.union([a,b],engine='manifold')
c.export(HERE/'diagnostic_mesh_union.stl')
np.savez_compressed(HERE/'diagnostic_mesh_union.npz',vertices=c.vertices,faces=c.faces)
(HERE/'diagnostic_mesh_union.json').write_text(json.dumps({'status':'Diagnostic mesh union ONLY; not native CAD fusion','watertight':c.is_watertight,'winding_consistent':c.is_winding_consistent,'volume_mm3':c.volume,'faces':len(c.faces),'elapsed_seconds':time.perf_counter()-started},indent=2)+'\n')
print((HERE/'diagnostic_mesh_union.json').read_text())
