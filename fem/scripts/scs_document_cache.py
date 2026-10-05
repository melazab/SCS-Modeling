"""Document-owned cache references and serialization-tolerant lead matching.

The numerical cache stays on disk. Reuse requires unchanged mesh parameters,
valid provenance manifests, and matching labelled triangle geometry. Vertex
numbering and float32 STL rounding are not physical geometry changes.
"""
import json
import os
from pathlib import Path
import numpy as np
import artifacts
import config as C
import live_lead
import stlio

GEOMETRY_TOL_MM = 2e-5  # 20 nm: covers binary STL rounding at this model's scale
_matches = {}


def same_params(a, b):
    def normalize(value):
        if isinstance(value, (float, int)):
            return round(value, 10)
        if isinstance(value, dict):
            return {k: normalize(v) for k, v in value.items()}
        if isinstance(value, (tuple, list)):
            return [normalize(v) for v in value]
        return value
    return normalize(a) == normalize(b)


def same_surface(current, saved):
    from scipy.spatial import cKDTree
    vertices, triangles = current
    other, faces = saved
    if vertices.shape != other.shape or triangles.shape != faces.shape:
        return False
    distances, mapping = cKDTree(other).query(vertices)
    if np.max(distances) > GEOMETRY_TOL_MM or len(np.unique(mapping)) != len(other):
        return False
    def canonical(t):
        t = np.sort(t, axis=1)
        return t[np.lexsort((t[:, 2], t[:, 1], t[:, 0]))]
    left = set(map(tuple, canonical(mapping[triangles])))
    right = set(map(tuple, canonical(faces)))
    if left == right:
        return True
    # OCCT can choose different diagonals on planar end caps after BRep
    # restore. Accept only coplanar patches with identical boundary edges.
    def planar_patches(tris):
        tris = list(tris)
        edges = {}
        for i, tri in enumerate(tris):
            for a,b in ((tri[0],tri[1]),(tri[1],tri[2]),(tri[0],tri[2])):
                edges.setdefault(tuple(sorted((a,b))), []).append(i)
        neighbors = [set() for _ in tris]
        xyz = other[np.array(tris)]
        normals = np.cross(xyz[:,1]-xyz[:,0], xyz[:,2]-xyz[:,0])
        lengths = np.linalg.norm(normals, axis=1)
        if np.any(lengths == 0):
            return None
        normals /= lengths[:,None]
        for owners in edges.values():
            if len(owners) > 2:
                return None
            for i in owners:
                for j in owners:
                    if np.max(np.abs((xyz[j]-xyz[i,0]) @ normals[i])) <= GEOMETRY_TOL_MM:
                        neighbors[i].add(j)
        unseen = set(range(len(tris)))
        patches = set()
        while unseen:
            todo = [unseen.pop()]; component = set(todo)
            while todo:
                for j in neighbors[todo.pop()] & unseen:
                    unseen.remove(j); component.add(j); todo.append(j)
            ids = np.unique([tris[i] for i in component])
            a,b,c = other[list(tris[next(iter(component))])]
            normal = np.cross(b-a,c-a); length = np.linalg.norm(normal)
            if length == 0 or np.max(np.abs((other[ids]-a) @ (normal/length))) > GEOMETRY_TOL_MM:
                return None
            boundary = tuple(sorted(edge for edge,owners in edges.items()
                                    if sum(i in component for i in owners)==1))
            patches.add(boundary)
        return patches
    lpatch = planar_patches(left-right)
    return lpatch is not None and lpatch == planar_patches(right-left)


def geometry_matches(arrays, run_dir):
    contacts, insulator = C.lead_stls(os.path.join(run_dir, 'lead'))
    mapping = getattr(arrays, 'contact_map', {})
    map_path = os.path.join(run_dir, 'lead', 'contact_map.json')
    if len({v[0] for v in mapping.values()}) > 1:
        try:
            with open(map_path) as f:
                if json.load(f) != {str(k): v for k, v in mapping.items()}:
                    return False
        except (OSError, ValueError):
            return False
    paths = dict(contacts, insulator=insulator)
    if set(paths) != set(arrays):
        return False
    try:
        key = (live_lead.fingerprint_from_arrays(arrays),
               tuple((str(k), artifacts.file_hash(p)) for k, p in paths.items()))
        if key not in _matches:
            _matches[key] = all(same_surface(arrays[k], stlio.read_stl(p)) for k,p in paths.items())
        return _matches[key]
    except (OSError, ValueError):
        return False


def run_params(run_dir):
    with open(os.path.join(run_dir, 'params.json')) as f:
        params = json.load(f)
    params.pop('LEAD_DIR', None)
    return params


def document_run(doc):
    group = doc.getObject('SCS_MeshPreviewGroup')
    if group and getattr(group, 'MeshRunDirectory', ''):
        return group.MeshRunDirectory
    montage = doc.getObject('SCS_Montage')
    if montage and getattr(montage, 'SolutionNpz', ''):
        return os.path.dirname(montage.SolutionNpz)
    return None


def document_params(doc):
    group = doc.getObject('SCS_MeshPreviewGroup')
    if group and getattr(group, 'MeshParameters', ''):
        return json.loads(group.MeshParameters)
    run = document_run(doc)
    if run:
        try:
            return run_params(run)
        except (OSError, ValueError):
            pass
    return None


def bind_mesh(doc, run_dir, params):
    group = doc.getObject('SCS_MeshPreviewGroup')
    if group is None:
        return
    for name in ('MeshRunDirectory', 'MeshParameters'):
        if name not in group.PropertiesList:
            group.addProperty('App::PropertyString', name, 'FEM cache')
        group.setEditorMode(name, 1)
    values = dict(MeshRunDirectory=os.path.abspath(run_dir), MeshParameters=json.dumps(params, sort_keys=True))
    for name,value in values.items():
        if getattr(group, name) != value:
            setattr(group, name, value)


def resolve_run(doc, arrays, params):
    direct = os.path.join(C.LEAD_RUNS, live_lead.run_key_for(live_lead.fingerprint_from_arrays(arrays), params))
    candidates = [document_run(doc), direct]
    # Migration for old mesh-only documents with no explicit cache reference.
    if not candidates[0]:
        candidates += [str(p.parent) for p in Path(C.LEAD_RUNS).glob('*/params.json')]
    for run in dict.fromkeys(candidates):
        if not run:
            continue
        try:
            if not same_params(run_params(run), params):
                continue
            if artifacts.preview_ready(run) and geometry_matches(arrays, run):
                bind_mesh(doc, run, params)
                return run
        except (OSError, ValueError, KeyError):
            continue
    return direct
