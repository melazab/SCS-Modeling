"""Mesh helpers for fusing and blending anatomical STL bodies.

Runs in src/extend/.venv (numpy, scipy, manifold3d, scikit-image), NOT in
fem/.venv: fem/requirements.txt is hashed by the solution provenance check.

    python3 -m venv src/extend/.venv
    src/extend/.venv/bin/pip install numpy scipy manifold3d scikit-image matplotlib

Launch anything that can grow under a memory cap, so only the job dies:
    systemd-run --user --scope -q -p MemoryMax=8G -p MemorySwapMax=0 <cmd>
"""
import os
import struct
import sys

import numpy as np
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "fem", "scripts"))
sys.path.insert(0, os.path.join(ROOT, "src", "ansys"))

import config as C  # noqa: E402  (read-only)
import stlio  # noqa: E402
from check_tissue_map import matches, parse_tissue_map  # noqa: E402
from inside import InsideTester  # noqa: E402,F401

ENTRIES = parse_tissue_map(C.TISSUE_MAP_YAML)
PIPE_OF_CLASS = {cls: pipe for pipe, cls in C.ANATOMY_CLASS}
ORDER = [p for p, _ in C.ANATOMY_CLASS]          # FEM priority, first wins


def pipe_of(stem):
    """fem/ tissue name of an STL stem (first tissue_map pattern wins), or None."""
    for t in ENTRIES:
        if any(matches(stem, p) for p in t["patterns"]):
            return PIPE_OF_CLASS.get(t["name"])
    return None


def read_tri(path):
    """(n, 3, 3) float64 triangles in file order."""
    with open(path, "rb") as f:
        raw = f.read()
    if raw[:5] == b"solid" and b"facet normal" in raw[:2048]:
        return stlio._read_ascii(raw)
    return stlio._read_binary(raw)


def write_stl(path, tri, header=b"src/extend"):
    tri = np.asarray(tri, dtype=np.float32)
    rec = np.zeros(len(tri), dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    ln = np.linalg.norm(n, axis=1, keepdims=True)
    rec["n"] = np.where(ln > 0, n / np.where(ln > 0, ln, 1), 0)
    rec["v"] = tri
    with open(path, "wb") as f:
        f.write(header[:80].ljust(80, b" "))
        f.write(struct.pack("<I", len(tri)))
        f.write(rec.tobytes())


def weld(tri, tol):
    """Indexed mesh; vertices closer than tol are merged (union-find over a
    k-d tree, so near-coincident pairs never straddle a rounding bin)."""
    P = tri.reshape(-1, 3)
    U, inv = np.unique(P, axis=0, return_inverse=True)
    inv = inv.ravel()
    if tol > 0:
        par = np.arange(len(U))

        def find(a):
            while par[a] != a:
                par[a] = par[par[a]]
                a = par[a]
            return a
        for a, b in cKDTree(U).query_pairs(tol):
            ra, rb = find(a), find(b)
            if ra != rb:
                par[max(ra, rb)] = min(ra, rb)
        root = np.array([find(a) for a in range(len(U))])
        keep, newi = np.unique(root, return_inverse=True)
        U, inv = U[keep], newi.ravel()[inv]
    T = inv.reshape(-1, 3)
    return U, T[(T[:, 0] != T[:, 1]) & (T[:, 1] != T[:, 2]) & (T[:, 0] != T[:, 2])]


def manifold(V, T):
    import manifold3d as m3
    return m3.Manifold(m3.Mesh(vert_properties=np.ascontiguousarray(V, np.float32),
                               tri_verts=np.ascontiguousarray(T, np.uint32)))


def load_manifold(path):
    """manifold3d solid of an STL, welding float cracks with the smallest
    tolerance that makes it a valid manifold. Returns (solid, tolerance)."""
    import manifold3d as m3
    tri = read_tri(path)
    for tol in (0, 1e-5, 1e-4, 1e-3):
        M = manifold(*weld(tri, tol))
        if M.status() == m3.Error.NoError:
            return M, tol
    raise RuntimeError("%s: not a closed manifold even welded at 1e-3 mm" % path)


def drop_slivers(M, eps=1e-6):
    """Remove zero-volume shells (and voids) left by booleans on coincident faces."""
    import manifold3d as m3
    keep = [c for c in M.decompose() if abs(c.volume()) > eps]
    return m3.Manifold.compose(keep) if keep else m3.Manifold()


def manifold_tri(M):
    mesh = M.to_mesh()
    V = np.asarray(mesh.vert_properties)[:, :3]
    return V[np.asarray(mesh.tri_verts)]


def tube_centreline(tri, step=0.5):
    """Axis of a tube from its side faces: each face centroid moved inward
    along its own (outward) normal by the radius, averaged per axial bin,
    area-weighted. Independent of tessellation; follows bends.
    Returns (points bottom->top, radius)."""
    P = tri.reshape(-1, 3)
    c = P.mean(0)
    ax = np.linalg.svd(P - c, full_matrices=False)[2][0]
    ax = ax if ax[2] > 0 else -ax
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    area = 0.5 * np.linalg.norm(n, axis=1)
    n = n / np.maximum(2 * area, 1e-30)[:, None]
    s_all = (P - c) @ ax
    mid = (s_all > s_all.min() + 2) & (s_all < s_all.max() - 2)
    r = float(np.median(np.linalg.norm((P - c)[mid] - np.outer(s_all[mid], ax), axis=1)))
    side = (np.abs(n @ ax) < 0.3) & (area > 0)
    q = tri.mean(1)[side] - r * n[side]
    sq, w = (q - c) @ ax, area[side]
    edges = np.arange(s_all.min(), s_all.max() + step, step)
    pts = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (sq >= a) & (sq < b)
        if w[m].sum() > 0:
            pts.append((q[m] * w[m, None]).sum(0) / w[m].sum())
    return np.array(pts), r


def tube_mesh(path, r, n=32):
    """Closed tube of radius r along a polyline (rotation-minimising frames,
    flat end caps), as (V, T) with outward orientation."""
    T_ = np.gradient(path, axis=0)
    T_ /= np.linalg.norm(T_, axis=1)[:, None]
    u = np.cross(T_[0], [1.0, 0, 0])
    if np.linalg.norm(u) < 1e-6:
        u = np.cross(T_[0], [0, 1.0, 0])
    U = [u / np.linalg.norm(u)]
    for i in range(1, len(path)):
        u = U[-1] - (U[-1] @ T_[i]) * T_[i]
        U.append(u / np.linalg.norm(u))
    U = np.array(U)
    W = np.cross(T_, U)
    th = 2 * np.pi * np.arange(n) / n
    ring = path[:, None, :] + r * (np.cos(th)[None, :, None] * U[:, None, :] + np.sin(th)[None, :, None] * W[:, None, :])
    m = len(path)
    V = np.vstack([ring.reshape(-1, 3), path[0], path[-1]])
    i, j = np.meshgrid(np.arange(m - 1), np.arange(n), indexing="ij")
    a, b = i * n + j, i * n + (j + 1) % n
    c_, d_ = (i + 1) * n + (j + 1) % n, (i + 1) * n + j
    F = np.vstack([np.column_stack([a.ravel(), b.ravel(), c_.ravel()]),
                   np.column_stack([a.ravel(), c_.ravel(), d_.ravel()])])
    c0, c1 = len(V) - 2, len(V) - 1
    jj = np.arange(n)
    F = np.vstack([F, np.column_stack([np.full(n, c0), (jj + 1) % n, jj]),
                   np.column_stack([np.full(n, c1), (m - 1) * n + jj, (m - 1) * n + (jj + 1) % n])])
    a_, b_, cc = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    if np.einsum("ij,ij->i", a_, np.cross(b_, cc)).sum() < 0:
        F = F[:, ::-1]
    return V, F
