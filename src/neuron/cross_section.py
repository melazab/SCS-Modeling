"""Exact planar cross-sections of STL bodies, as closed polygons.

Panel A of the poster figure used to be a rasterised point-in-shell mask on a
0.15 mm grid, so every boundary was a pixel staircase and the lead had to be
faked with rectangles. This computes the real thing: intersect each triangle
with the plane z = z0, then chain the resulting segments into closed loops.
The output is vector geometry -- crisp at any print size, and it renders the
lead's true (elliptical, because the lead is tilted) profile instead of a box.

Degenerate cases are avoided by nudging the plane off any vertex it lands on,
rather than by special-casing vertex-on-plane triangles.
"""
import numpy as np


def _segments(v, t, z):
    """Plane/triangle intersection segments at height z, as (m, 2, 2) xy pairs."""
    tri = v[t]                                   # (n, 3, 3)
    d = tri[:, :, 2] - z                         # signed height of each vertex
    above = d > 0
    n_above = above.sum(axis=1)
    cross = (n_above == 1) | (n_above == 2)      # triangle straddles the plane
    if not cross.any():
        return np.empty((0, 2, 2))
    tri, d, above = tri[cross], d[cross], above[cross]
    pts = []
    for e0, e1 in ((0, 1), (1, 2), (2, 0)):      # the three edges
        sel = above[:, e0] != above[:, e1]       # this edge crosses the plane
        a, b = tri[:, e0], tri[:, e1]
        da, db = d[:, e0], d[:, e1]
        w = np.where(sel, da / np.where(da - db == 0, np.nan, da - db), np.nan)
        p = a[:, :2] + (b[:, :2] - a[:, :2]) * w[:, None]
        pts.append(np.where(sel[:, None], p, np.nan))
    P = np.stack(pts, axis=1)                    # (m, 3, 2), one row is NaN
    good = np.isfinite(P).all(axis=2)
    out = np.empty((len(P), 2, 2))
    for i in range(len(P)):
        out[i] = P[i][good[i]][:2]
    return out


def _chain(seg, tol=2e-4):
    """Join segments end-to-end into closed loops.

    Endpoints are welded with a KD-tree rather than matched on rounded keys:
    the same physical point reached from two adjacent triangles can differ by
    ~1e-5 mm in an STL's float coordinates, which a tight key silently splits,
    leaving loops open and holes unfilled.
    """
    if len(seg) == 0:
        return []
    from scipy.spatial import cKDTree
    pts = seg.reshape(-1, 2)
    tree = cKDTree(pts)
    # Weld: every point adopts the lowest index within tol of it.
    node = np.arange(len(pts))
    for a, b in tree.query_pairs(tol):
        lo, hi = (a, b) if node[a] < node[b] else (b, a)
        node[node == node[hi]] = node[lo]
    e = node.reshape(-1, 2)
    adj = {}
    for k, (a, b) in enumerate(e):
        if a != b:
            adj.setdefault(a, []).append((b, k))
            adj.setdefault(b, []).append((a, k))
    used, loops = set(), []
    for seed in adj:
        for _, k0 in adj[seed]:
            if k0 in used:
                continue
            loop, cur, edge = [pts[2 * k0]], seed, k0
            while True:
                used.add(edge)
                nxt = e[edge][1] if e[edge][0] == cur else e[edge][0]
                loop.append(pts[2 * edge + (1 if e[edge][0] == cur else 0)])
                if nxt == seed:
                    break
                step = next(((n, kk) for n, kk in adj.get(nxt, ()) if kk not in used), None)
                if step is None:
                    break
                cur, edge = nxt, step[1]
            if len(loop) >= 4:
                loops.append(np.array(loop))
    return loops


def section(v, t, z, snap=1e-4):
    """Closed xy loops where the body (v, t) meets the plane z.

    The plane is nudged off any vertex that lies exactly on it, which removes
    the coplanar-triangle and vertex-touch degeneracies outright.
    """
    zs = v[:, 2]
    while np.any(np.abs(zs - z) < snap):
        z += snap
    return _chain(_segments(v, t, z))
