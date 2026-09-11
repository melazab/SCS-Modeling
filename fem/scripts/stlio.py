"""Minimal binary/ASCII STL reader + topology audit, numpy only.

Used by the FreeCAD-side FEM pipeline (fem/). Reads the RADO STL shells from
STL_files/ -- every body in NBF_RADO-SCS_fem.FCStd has an identity Placement,
so STL coordinates ARE document coordinates and the generated lead STLs in
fem/leads/ drop straight in untransformed.
"""
import numpy as np
import struct


def read_stl(path):
    """Return (verts (n,3) float64 mm, tris (m,3) int) with vertices welded."""
    with open(path, "rb") as f:
        head = f.read(5)
        f.seek(0)
        raw = f.read()
    if head == b"solid" and b"facet normal" in raw[:2048]:
        tri = _read_ascii(raw)
    else:
        tri = _read_binary(raw)
    return weld(tri)


def _read_binary(raw):
    n = struct.unpack("<I", raw[80:84])[0]
    body = np.frombuffer(raw[84:84 + n * 50], dtype=np.uint8)
    body = body.reshape(n, 50)
    f = body[:, 12:48].copy().view("<f4").reshape(n, 3, 3)
    return f.astype(np.float64)


def _read_ascii(raw):
    vals = []
    for line in raw.decode("ascii", "replace").splitlines():
        s = line.strip()
        if s.startswith("vertex"):
            vals.append([float(x) for x in s.split()[1:4]])
    a = np.asarray(vals, dtype=np.float64)
    return a.reshape(-1, 3, 3)


def weld(tri, tol=1e-7):
    """Weld coincident vertices onto a shared index array."""
    pts = tri.reshape(-1, 3)
    q = np.round(pts / tol).astype(np.int64)
    _, first, inv = np.unique(q, axis=0, return_index=True, return_inverse=True)
    verts = pts[first]
    tris = inv.reshape(-1, 3)
    # drop degenerate triangles (two welded corners identical)
    good = (tris[:, 0] != tris[:, 1]) & (tris[:, 1] != tris[:, 2]) & (tris[:, 0] != tris[:, 2])
    return verts, tris[good]


def audit(verts, tris):
    """Topology report: edge manifoldness, closure, shell count, volume."""
    e = np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]])
    key = np.sort(e, axis=1)
    uniq, counts = np.unique(key, axis=0, return_counts=True)
    # signed volume via divergence theorem
    a, b, c = verts[tris[:, 0]], verts[tris[:, 1]], verts[tris[:, 2]]
    vol = np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0
    area = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
    return dict(
        nv=len(verts), nt=len(tris),
        boundary_edges=int((counts == 1).sum()),
        nonmanifold_edges=int((counts > 2).sum()),
        closed=bool((counts == 2).all()),
        volume_mm3=float(vol),
        n_shells=n_shells(tris, len(verts)),
        min_area=float(area.min()), degenerate=int((area < 1e-9).sum()),
        bbox=(verts.min(0).tolist(), verts.max(0).tolist()),
    )


def n_shells(tris, nv):
    """Connected components over the vertex graph (union-find)."""
    parent = np.arange(nv)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for t in tris:
        r0 = find(t[0])
        for v in t[1:]:
            r1 = find(v)
            if r0 != r1:
                parent[r1] = r0
    used = np.unique(tris)
    return len({find(v) for v in used})
