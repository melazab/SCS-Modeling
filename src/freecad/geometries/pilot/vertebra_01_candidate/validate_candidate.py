"""Reproducible bidirectional area sampling; estimates, not Hausdorff proof."""
from pathlib import Path
import json,time,sys
import numpy as np
import trimesh
HERE=Path(__file__).resolve().parent; BASE=HERE.parent/'vertebra_01'
start=time.perf_counter(); np.random.seed(20261005)
data=json.loads((BASE/'measurements.json').read_text()); R=np.array(data['rotation_global_to_local'])
source=trimesh.load_mesh(BASE/'seam_closed_reference.stl')
filename=sys.argv[1] if len(sys.argv)>1 else 'body_candidate.stl'
is_hybrid=filename!='body_candidate.stl'
if filename.endswith('.npz'):
 raw=np.load(HERE/filename); candidate=trimesh.Trimesh(raw['vertices'],raw['faces'],process=False)
else: candidate=trimesh.load_mesh(HERE/filename)
def distances(target,p):
 vals=[]
 for i in range(0,len(p),1000):
  vals.extend(trimesh.proximity.closest_point(target,p[i:i+1000])[1])
 return np.asarray(vals)
res={'mode':filename+' whole surface' if is_hybrid else 'body region: local y < 39 mm only','samples_requested_per_direction':100000,'candidate_mesh_watertight':candidate.is_watertight,'candidate_mesh_winding_consistent':candidate.is_winding_consistent,'candidate_mesh_volume_mm3':candidate.volume,'source_closed_reference_volume_mm3':source.volume,'relative_volume_error':candidate.volume/source.volume-1 if is_hybrid else None,'candidate_mesh_faces':len(candidate.faces)}
all_d=[]
for name,a,b in [('source_to_candidate',source,candidate),('candidate_to_source',candidate,source)]:
 p,_=trimesh.sample.sample_surface(a,100000,seed=20261005)
 if not is_hybrid: p=p[(p@R.T)[:,1]<39]
 d=distances(b,p); all_d.extend(d)
 res[name]={'n_area_samples':len(p),'rms_mm':float(np.sqrt(np.mean(d*d))),'sample_max_mm':float(d.max()),'p95_mm':float(np.percentile(d,95)),'max_point_mm':p[d.argmax()].tolist()}
 print(name,res[name],flush=True)
res['pooled_rms_mm']=float(np.sqrt(np.mean(np.array(all_d)**2)))
res['bidirectional_sample_max_mm']=float(max(all_d))
res['runtime_seconds']=time.perf_counter()-start
(HERE/('validation_union.json' if is_hybrid else 'validation_body.json')).write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(res,indent=2))
