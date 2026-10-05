"""Content provenance and atomic publication for mesh/solution artifacts.

Unmarked files are historical results, never automatic cache hits. File hashes
are memoized by inode/size/mtime/ctime so GUI polling does not reread all STLs.
"""
import hashlib
import ast
import json
import os
from pathlib import Path
import tempfile

import config as C

_HASHES = {}
_CLASSIFICATION_HASH = None


def file_hash(path):
    path = os.path.abspath(path)
    s = os.stat(path)
    stamp = (s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)
    if path not in _HASHES or _HASHES[path][0] != stamp:
        h = hashlib.sha256()
        with open(path, 'rb') as f:
            for block in iter(lambda: f.read(1024 * 1024), b''):
                h.update(block)
        _HASHES[path] = (stamp, h.hexdigest())
    return _HASHES[path][1]


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def model_signature(stage='mesh'):
    global _CLASSIFICATION_HASH
    paths = [C.TISSUE_MAP_YAML, os.path.join(C.ROOT, 'fem', 'requirements.txt')]
    paths += [p for group in C.TISSUE_BODIES.values() for p in group]
    here = Path(__file__).parent
    names = ['artifacts.py', 'config.py', 'stlio.py', 'inside.py',
             'build_mesh.py', 'mesh_geometry.py', 'mesh_faces.py', 'mesh_preview.py']
    if stage == 'solve':
        names += ['assign_and_solve.py', 'solve_lead.py', 'element_batches.py',
                  'elmer_backend.py', 'solve_resources.py']
        paths += [p for p in ('/opt/elmerfem/bin/ElmerSolver_mpi',
                             '/opt/elmerfem/bin/ElmerGrid',
                             '/opt/elmerfem/lib/elmersolver/libelmersolver.so',
                             '/opt/elmerfem/share/elmersolver/lib/StatCurrentSolve.so',
                             '/opt/elmerfem/deps/usr/lib/x86_64-linux-gnu/libHYPRE.so')
                  if os.path.isfile(p)]
    paths += [str(here / name) for name in names]
    # Classification shares a module with the solver. Hash its actual code
    # separately, so changing CG/AMG does not invalidate a completed mesh.
    assignment = here / 'assign_and_solve.py'
    stamp = file_hash(assignment)
    if _CLASSIFICATION_HASH is None or _CLASSIFICATION_HASH[0] != stamp:
        tree = ast.parse(assignment.read_text())
        selected = [ast.dump(node, include_attributes=False) for node in tree.body
                    if isinstance(node, ast.FunctionDef)
                    and node.name in ('build_testers', 'classify_tets', 'dura_leak_report')]
        _CLASSIFICATION_HASH = (stamp, digest(selected))
    return digest(dict(files=[(os.path.relpath(p, C.ROOT), file_hash(p)) for p in sorted(paths)],
                       classification=_CLASSIFICATION_HASH[1]))


def mesh_signature(params, contact_stl, insulator_stl):
    return digest(dict(model=model_signature(), params=params,
                       lead=[(str(i), file_hash(p)) for i, p in sorted(contact_stl.items())]
                       + [('insulator', file_hash(insulator_stl))]))


def solution_signature(run_dir, background=False):
    return digest(dict(model=model_signature('solve'), mesh=file_hash(os.path.join(run_dir, 'mesh.npz')),
                       background=bool(background)))


def atomic_json(path, value):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(os.path.abspath(path)), suffix='.tmp')
    try:
        with os.fdopen(fd, 'w') as f:
            json.dump(value, f, indent=2)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def atomic_npz(path, **arrays):
    import numpy as np
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(os.path.abspath(path)), suffix='.tmp')
    try:
        with os.fdopen(fd, 'wb') as f:
            np.savez_compressed(f, **arrays)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def publish(path, signature):
    atomic_json(path + '.manifest.json', dict(signature=signature, sha256=file_hash(path)))


def valid(path, signature):
    try:
        with open(path + '.manifest.json') as f:
            m = json.load(f)
        return m['signature'] == signature and m['sha256'] == file_hash(path)
    except (OSError, ValueError, KeyError):
        return False


def preview_ready(run_dir):
    try:
        with open(os.path.join(run_dir, 'mesh_report.json')) as f:
            report = json.load(f)
        sig = report['signature']
        return (report['model_signature'] == model_signature()
                and valid(os.path.join(run_dir, 'mesh.npz'), sig)
                and valid(os.path.join(run_dir, 'preview.vtp'), sig))
    except (OSError, ValueError, KeyError):
        return False


def result_path(name):
    """Keep the background scenario's derived artifacts separate."""
    if os.environ.get('SCS_BACKGROUND', '0') == '1':
        stem, ext = os.path.splitext(name)
        name = stem + '_bg' + ext
    return os.path.join(C.OUT, name)
