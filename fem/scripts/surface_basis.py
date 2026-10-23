"""Per-contact potentials at every surface the Potential Visualizer draws, so a
cluster solve can stay on the cluster.

    python surface_basis.py <run_dir> [--background]

A solution.npz for the 415 M-tet T6-T12 mesh is ~12 GB (the whole mesh plus
one field per contact), hours over the VPN. Plotting only ever reads the field
at the vertices of the tissue and lead surfaces (scs_montage.sample_tissue),
and any montage is a weighted sum of the per-contact fields. So this samples
each contact's field at those vertices once, next to the solution, and
writes surface_basis.npz (surface_basis_bg.npz for solution_bg.npz):

    names            the surfaces, as scs_montage.surface_paths() keys them
    phi_<name>       (contacts, vertices) V/A
    E_<name>         (contacts, vertices, 3) V/m per A, float32, from the
                     tissue's OWN tetrahedra (lead surfaces: the tissue next
                     to them). E's normal component jumps at an interface by
                     the conductivity ratio (CSF 1.7 vs white matter 0.083
                     S/m across the fibres), so a vertex taken from the
                     neighbour's side could read ~20x low.
    verts_<name>     sha256 of the vertex array it was sampled at
    filled_<name>    vertices taken from the nearest solved node
    contact_ids, solution_sha256, self_check_ok, self_check_msg

Sampling is sample_tissue()'s, contact by contact: barycentric inside the
located tetrahedron, and for a vertex the strict containment test misses, the
nearest solved node's value. A montage from this file equals the one from the
solution itself (scs_montage reads it when the solution is not here). The
self check needs the whole solution, so it runs here and its verdict is kept.
"""
import argparse
import hashlib
import os
import sys
import time

import numpy as np

import artifacts
import scs_montage as M

HERE = os.path.dirname(os.path.abspath(__file__))


def basis_path(sol_npz):
    """surface_basis.npz beside solution.npz, surface_basis_bg.npz beside solution_bg.npz."""
    head, name = os.path.split(sol_npz)
    return os.path.join(head, 'surface_basis_bg.npz' if name == 'solution_bg.npz' else 'surface_basis.npz')


def verts_hash(verts):
    return hashlib.sha256(np.ascontiguousarray(verts, dtype='<f8').tobytes()).hexdigest()


def signature(solution_signature):
    """Ties the basis to its solution and to the code that sampled it."""
    return artifacts.digest(dict(solution=solution_signature,
                                 code=[artifacts.file_hash(os.path.join(HERE, n))
                                       for n in ('surface_basis.py', 'field.py')]))


def valid(sol_npz):
    """The basis beside sol_npz matches that solution (by its manifest, so the
    solution itself need not be here) and this code."""
    import live_lead
    run_dir = os.path.dirname(sol_npz)
    try:
        sol_sig = live_lead.solution_signature(run_dir, os.path.basename(sol_npz) == 'solution_bg.npz')
    except (OSError, ValueError, KeyError):
        return False
    return artifacts.valid(basis_path(sol_npz), signature(sol_sig))


def tissue_E(f, keep, verts):
    """Per-contact E (V/m per A; coordinates are mm) at verts from the
    tetrahedra where keep is True: the one containing each vertex, else the
    one whose centroid is nearest. A P1 field has one gradient per
    tetrahedron: grad = A^-1 dphi, A's rows the edges from node 0."""
    from field import TetField
    if not keep.any():                                      # tissue not in this (cropped) mesh
        return np.full((f.values.shape[0], len(verts), 3), np.nan, np.float32), 0
    sub = TetField(f.nodes, f.tets[keep], f.values)
    ti, _w = sub.locate(verts)
    miss = ti < 0
    if miss.any():
        ti[miss] = sub.tree.query(verts[miss])[1]
    tt = sub.tets[ti]
    p = sub.nodes[tt]
    A = p[:, 1:] - p[:, :1]                                 # (n, 3, 3)
    dphi = (sub.values[:, tt[:, 1:]] - sub.values[:, tt[:, :1]]).transpose(1, 2, 0)   # (n, 3, contacts)
    grad = np.linalg.solve(A, dphi)                         # (n, 3, contacts) V/mm per A
    return (-1000.0 * grad).transpose(2, 0, 1).astype(np.float32), int(miss.sum())


def build(run_dir, background=False):
    import live_lead
    sol = os.path.join(run_dir, 'solution_bg.npz' if background else 'solution.npz')
    sol_sig = live_lead.solution_signature(run_dir, background)
    if not artifacts.valid(sol, sol_sig):
        raise SystemExit('%s does not match its manifest' % sol)
    t0 = time.time()
    f, d, tree = M._load_basis(sol)
    print('loaded %s: %d contacts, %d nodes (%.0f s)' % (os.path.basename(sol), *f.values.shape,
                                                         time.time() - t0), flush=True)
    out = dict(contact_ids=d['contact_ids'], solution_sha256=artifacts.file_hash(sol))
    lab, order = d['label'], [str(n) for n in d['order']]
    lead = np.isin(lab, [i for i, n in enumerate(order) if n.startswith(('contact', 'insulator'))])
    names = []
    for name, path in M.surface_paths(sol).items():
        verts, _tris = M.read_surface(path)
        phi = f(verts)                                   # (contacts, vertices), NaN outside
        miss = np.isnan(phi).any(axis=0)
        if miss.any():
            dist, idx = tree.query(verts[miss])
            phi[:, miss] = f.values[:, idx]
            far = int((dist > M.FILL_TOL_MM).sum())
            print('  %s: %d of %d vertices from the nearest node (max %.3f mm%s)'
                  % (name, miss.sum(), len(verts), dist.max(),
                     '' if not far else '; %d beyond %.1f mm' % (far, M.FILL_TOL_MM)), flush=True)
        out['phi_' + name] = phi
        own = name if name in order else name.replace('_', '')     # contact_3 -> contact3
        keep = (lab == order.index(own)) if own in order and not own.startswith(('contact', 'insulator')) \
            else (~lead & (lab >= 0))
        out['E_' + name], far = tissue_E(f, keep, verts)
        if far:
            print('  %s: E for %d vertices from the nearest tetrahedron of its tissue' % (name, far), flush=True)
        out['verts_' + name] = verts_hash(verts)
        out['filled_' + name] = int(miss.sum())
        names.append(name)
    ok, msg = M.self_check(sol_npz=sol)
    print('self check: %s, %s' % ('ok' if ok else 'FAILED', msg), flush=True)
    out.update(names=np.array(names), self_check_ok=ok, self_check_msg=msg)
    path = basis_path(sol)
    artifacts.atomic_npz(path, **out)
    artifacts.publish(path, signature(sol_sig))
    print('wrote %s: %d surfaces, %.0f MB (%.0f s)' % (path, len(names), os.path.getsize(path)/1e6,
                                                      time.time() - t0), flush=True)
    return path


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('run_dir')
    p.add_argument('--background', action='store_true')
    a = p.parse_args()
    build(os.path.abspath(a.run_dir), a.background)
    sys.exit(0)
