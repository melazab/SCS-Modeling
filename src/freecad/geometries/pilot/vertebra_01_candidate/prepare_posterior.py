"""Preserve double precision when passing a cut mesh into the CAD kernel."""
import json
from pathlib import Path
import numpy as np
import trimesh
HERE=Path(__file__).resolve().parent; BASE=HERE.parent/'vertebra_01'
data=json.loads((BASE/'measurements.json').read_text())
rotation=np.array(data['rotation_global_to_local'])
mesh=trimesh.load_mesh(BASE/'seam_closed_reference.stl')
mesh.vertices=mesh.vertices@rotation.T
mesh=mesh.slice_plane([0,39.5,0],[0,1,0],cap=True)
mesh.vertices=mesh.vertices@rotation
assert mesh.is_watertight and mesh.is_winding_consistent
(HERE/'posterior_working.json').write_text(json.dumps({'vertices':mesh.vertices.tolist(),'faces':mesh.faces.tolist()}))
print({'faces':len(mesh.faces),'minimum_triangle_area_mm2':float(mesh.area_faces.min()),'volume_mm3':mesh.volume})
