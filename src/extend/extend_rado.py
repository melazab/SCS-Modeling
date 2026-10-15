"""Extend the RADO-SCS geometry rostrally by whole vertebral levels.

    fem/.venv/bin/python src/extend/extend_rado.py              # 4 levels: T6-T9
    fem/.venv/bin/python src/extend/extend_rado.py --levels 3

Writes fem/out/geometry/rado_up<N>/ (git-ignored): a complete, self-contained
STL set (RADO's files, the new levels, and the few RADO bodies this replaces),
manifest.json (inputs with SHA-256, every parameter, transform and check) and
extension_check.png. STL_files/ is only read.

WHY THIS IS POSSIBLE (measured, 2026-10-09). RADO is one level unit, copied:
224 of its 232 per-level bodies are exact rigid copies of a sibling at another
level (residual ~2e-5 mm), and the 8 bodies that run the whole length (white,
grey, CSF, meninges, epidural, two blood layers, aorta) are constant
cross-sections swept along a curve. RADO's levels were each placed by hand,
so the new levels copy RADO's TOP level arrangement (top vertebra, the disc
above it, the roots/DRG/sleeves/vessels/ganglia at that disc), not an average.

HOW THE LEVELS ARE STACKED. The step from RADO's middle to its top vertebra
(V1-3 -> V1-4, T11 -> T10 in disc naming) is a rigid motion S: a 5.87 deg
rotation about a left-right axis ~31 cm anterior to the vertebral bodies.
Level k is S^k applied to the top unit. The vertebrae therefore continue
RADO's own curve, a circular arc through T11 and T10, and each new disc sits
on and under its vertebrae exactly as RADO's top disc sits on T10. (RADO's
lower step, T12 -> T11, is 4.61 deg; the curve has no other definition, no
spline.) An earlier version continued straight and gave the anterior
vertebral edge a sawtooth; Mohamed's review rejected it, 2026-10-12.

The continuous bodies follow the same motion: the section at the join is
carried by S^t for real t (the S-orbit), so a copied vertebra's canal,
foramina and root entries meet the swept canal contents exactly as at T10.
RADO's own cord is not quite S-invariant: at the cut its axis is 1.5 deg off
the orbit, so the extension turns by 1.5 deg there.

ONE BODY PER CONTINUOUS TISSUE. Each continuous body is written as one closed
surface: the original cut at one plane shared by all eight, walls swept from
the cut points themselves, and a top cap. The cut and the walls share vertices, so the result
is closed and consistently oriented; that is checked, not assumed. It replaces
the original file in the set.

THE SYMPATHETIC TRUNK. RADO's trunk segments (neuro_sympathetic_a, ~29 mm
tubes, r 0.75 mm) were placed by hand and overlap their neighbours by
0.7-1.1 mm. They are not S-periodic: S leaves a copied segment 1.1 mm (right)
and 1.8 mm (left) off the one below. So each segment's upper end is bent, with
a smoothstep from its ganglion up to its top end, until it overlaps the
segment above by CHAIN_OVERLAP_MM. Ganglia and rami stay rigid. RADO's own top
segments get the same bend (written as <stem>__joined, replacing the
original), and the topmost copy is left free.

OVERLAPS. RADO's parts were never cut against each other, so its bodies
overlap (18 cm^3 between different tissues; src/extend/overlap_audit.py). The
copies carry the same overlaps, and fem/ resolves them by tissue priority, as
it always has for RADO. Removing them is left to the parametric model.

Level names follow RADO's disc filenames (top vertebra = T10, so 4 new levels
are T9..T6), the same convention as src/freecad/README.md "Vertebral levels";
the manifest also records the paper's naming (one level higher). New discs are
named for the levels they separate (T_disk_8_9x...), other copies keep their
source name plus "__up<k>"; every new name is asserted to map to the same
tissue_map.yaml class as its source.
"""
import argparse
import hashlib
import json
import os
import re
import shutil
import struct
import sys
import time

import numpy as np
from scipy.linalg import expm, logm
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "fem", "scripts"))
sys.path.insert(0, os.path.join(ROOT, "src", "ansys"))

import config as C  # noqa: E402  (read-only: hashed by the provenance check)
from check_tissue_map import matches, parse_tissue_map  # noqa: E402
from inside import InsideTester  # noqa: E402
import stlio  # noqa: E402

PIPE_OF_CLASS = {cls: pipe for pipe, cls in C.ANATOMY_CLASS}
CONTINUOUS_SPAN_MM = 80.0      # z-span above which a body is "continuous"
RIGID_TOL_MM = 1e-3
CHAIN_OVERLAP_MM = 0.8         # RADO's own trunk segments overlap 0.74-1.09 mm
SWEEP_MAX_DEG = 0.2            # rotation between consecutive swept rings


# --------------------------------------------------------------------------
# STL I/O and helpers
# --------------------------------------------------------------------------
def read_raw(path):
    """(n, 3, 3) float64 triangles in FILE order (welding reorders vertices,
    which hides exact rigid copies; raw order preserves them)."""
    with open(path, "rb") as f:
        raw = f.read()
    if raw[:5] == b"solid" and b"facet normal" in raw[:2048]:
        return stlio._read_ascii(raw)
    return stlio._read_binary(raw)


def write_stl(path, tri, header=b"src/extend/extend_rado.py"):
    tri = np.asarray(tri, dtype=np.float64)
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    ln = np.linalg.norm(n, axis=1, keepdims=True)
    n = np.where(ln > 0, n / np.where(ln > 0, ln, 1), 0)
    rec = np.zeros(len(tri), dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
    rec["n"], rec["v"] = n, tri
    with open(path, "wb") as f:
        f.write(header[:80].ljust(80, b" "))
        f.write(struct.pack("<I", len(tri)))
        f.write(rec.tobytes())


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def class_of(stem, entries):
    """tissue_map.yaml class, first pattern match wins (as config.py does)."""
    for t in entries:
        if any(matches(stem, p) for p in t["patterns"]):
            return t["name"]
    return None


def rigid_fit(A, B):
    """Proper rigid (R, t) minimising |R A + t - B|, and the max residual."""
    ca, cb = A.mean(0), B.mean(0)
    U, _, Vt = np.linalg.svd((A - ca).T @ (B - cb))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    t = cb - R @ ca
    return R, t, float(np.linalg.norm(A @ R.T + t - B, axis=1).max())


def apply(M, pts):
    return pts @ M[:3, :3].T + M[:3, 3]


def topology(V, T):
    """Closure and orientation of an indexed surface: boundary and
    non-manifold edge counts, and edges whose two faces disagree on direction
    (a consistently oriented closed surface has none)."""
    e = np.concatenate([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]])
    key = np.sort(e, axis=1)
    _, inv, counts = np.unique(key, axis=0, return_inverse=True, return_counts=True)
    fwd = np.bincount(inv.ravel(), weights=(e[:, 0] < e[:, 1]), minlength=len(counts))
    two = counts == 2
    a, b, c = V[T[:, 0]], V[T[:, 1]], V[T[:, 2]]
    return {"triangles": int(len(T)), "boundary_edges": int((counts == 1).sum()),
            "nonmanifold_edges": int((counts > 2).sum()),
            "misoriented_edges": int((two & (fwd != 1)).sum()),
            "volume_mm3": round(float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0), 3)}


class Step:
    """A rigid motion S (4x4) and its real powers S^t = expm(t log S): the
    screw motion that interpolates it. RADO's vertebral step is a rotation
    about a left-right axis (its translation along that axis is 0.04 mm)."""

    def __init__(self, M):
        self.M = M
        self.L = np.real(logm(M))
        assert np.abs(expm(self.L) - M).max() < 1e-9

    def power(self, t):
        if float(t).is_integer():
            return np.linalg.matrix_power(self.M, int(t)) if t >= 0 else \
                np.linalg.matrix_power(np.linalg.inv(self.M), int(-t))
        return expm(t * self.L)

    def angle_deg(self):
        return float(np.degrees(np.linalg.norm([self.L[2, 1], self.L[0, 2], self.L[1, 0]])))

    def axis(self):
        w = np.array([self.L[2, 1], self.L[0, 2], self.L[1, 0]])
        return w / np.linalg.norm(w)

    def speed(self, p):
        """mm moved per unit t by point p (= arc length per level at p)."""
        return float(np.linalg.norm(self.L[:3, :3] @ p + self.L[:3, 3]))

    def centre_yz(self):
        R, t = self.M[:3, :3], self.M[:3, 3]
        return np.linalg.solve(np.eye(2) - R[1:, 1:], t[1:])


# --------------------------------------------------------------------------
# Centreline
# --------------------------------------------------------------------------
def cord_centroids(paths, step=1.0, h=0.1):
    """z-slice centroids of white+grey matter (an oblique planar cut of a
    tube has its centroid on the axis, so tilt does not bias this)."""
    testers, verts = [], []
    for p in paths:
        v, t = stlio.read_stl(p)
        testers.append(InsideTester(v, t))
        verts.append(v)
    V = np.vstack(verts)
    xs = np.arange(V[:, 0].min() - 0.5, V[:, 0].max() + 0.5, h)
    zs = np.arange(V[:, 2].min() + 3.0, V[:, 2].max() - 3.0, step)
    out, yc = [], None
    for z in zs:
        ylo, yhi = (V[:, 1].min() - .5, V[:, 1].max() + .5) if yc is None else (yc - 6, yc + 6)
        ys = np.arange(ylo, yhi, h)
        X, Y = np.meshgrid(xs, ys, indexing="ij")
        pts = np.column_stack([X.ravel(), Y.ravel(), np.full(X.size, z)])
        m = np.zeros(len(pts), bool)
        for t in testers:
            m |= t(pts)
        c = pts[m].mean(0)
        yc = c[1]
        out.append(c)
    return np.array(out)


class Curve:
    """Planar centreline in the y-z plane at x = x0, as dense samples of arc
    length s, (y, z) and tilt phi (tangent = (0, sin phi, cos phi))."""

    def __init__(self, x0, Y, Z):
        self.x0 = x0
        self.Y, self.Z = np.asarray(Y, float), np.asarray(Z, float)
        self.S = np.concatenate([[0], np.cumsum(np.hypot(np.diff(self.Y), np.diff(self.Z)))])
        self.PHI = np.arctan2(np.gradient(self.Y), np.gradient(self.Z))
        self._tree = cKDTree(np.column_stack([self.Y, self.Z]))

    @classmethod
    def cubic(cls, coef, x0, z_lo, z_hi, ds=0.05):
        z = np.arange(z_lo, z_hi + 1e-9, ds)
        return cls(x0, np.polyval(coef, z), z)

    def continued(self, s_join, step, t_end, ds=0.05):
        """This curve up to s_join, then the S-orbit of its point there."""
        keep = self.S < s_join
        p0 = self.point(s_join)
        ts = np.linspace(0, t_end, max(2, int(np.ceil(t_end * step.speed(p0) / ds))))
        P = np.array([apply(step.power(t), p0[None])[0] for t in ts])
        return Curve(self.x0, np.concatenate([self.Y[keep], P[:, 1]]), np.concatenate([self.Z[keep], P[:, 2]]))

    def point(self, s):
        return np.array([self.x0, np.interp(s, self.S, self.Y), np.interp(s, self.S, self.Z)])

    def phi(self, s):
        return float(np.interp(s, self.S, self.PHI))

    def tangent(self, s):
        p = self.phi(s)
        return np.array([0.0, np.sin(p), np.cos(p)])

    def normal(self, s):
        """In-plane unit normal, pointing toward +y (posterior)."""
        p = self.phi(s)
        return np.array([0.0, np.cos(p), -np.sin(p)])

    def project(self, pts):
        """Arc length of the nearest centreline sample, by (y, z)."""
        pts = np.atleast_2d(pts)
        _, i = self._tree.query(pts[:, 1:3])
        return self.S[i]


# --------------------------------------------------------------------------
# Cross-sections, cutting, sweeping, capping
# --------------------------------------------------------------------------
def slice_loops(V, T, p0, nrm):
    """Closed polylines where the plane (p0, nrm) cuts the welded mesh."""
    d = (V - p0) @ nrm
    d = np.where(np.abs(d) < 1e-9, 1e-9, d)
    sg = d > 0
    cross = np.flatnonzero(~(sg[T].all(1) | (~sg[T]).all(1)))
    edge_pt, seg = {}, []
    for ti in cross:
        a, b, c = T[ti]
        ends = []
        for u, w in ((a, b), (b, c), (c, a)):
            if sg[u] != sg[w]:
                key = (min(u, w), max(u, w))
                if key not in edge_pt:
                    edge_pt[key] = V[u] + d[u] / (d[u] - d[w]) * (V[w] - V[u])
                ends.append(key)
        if len(ends) == 2:
            seg.append(ends)
    nbr = {}
    for e0, e1 in seg:
        nbr.setdefault(e0, []).append(e1)
        nbr.setdefault(e1, []).append(e0)
    seen, loops = set(), []
    for start in nbr:
        if start in seen:
            continue
        loop, prev, cur = [start], None, start
        seen.add(start)
        while True:
            nx = [e for e in nbr[cur] if e != prev]
            if not nx:
                break
            nxt = nx[0]
            if nxt == start:
                break
            if nxt in seen:
                break
            loop.append(nxt)
            seen.add(nxt)
            prev, cur = cur, nxt
        loops.append(np.array([edge_pt[k] for k in loop]))
    return loops


def clean_loop(P, tol=1e-6):
    """Drop repeated and collinear points of a closed 2D polyline."""
    keep = [P[0]]
    for p in P[1:]:
        if np.linalg.norm(p - keep[-1]) > tol:
            keep.append(p)
    if len(keep) > 1 and np.linalg.norm(keep[0] - keep[-1]) <= tol:
        keep.pop()
    P = np.array(keep)
    changed = True
    while changed and len(P) > 3:
        a, b, c = np.roll(P, 1, 0), P, np.roll(P, -1, 0)
        cr = np.abs((b[:, 0] - a[:, 0]) * (c[:, 1] - b[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - b[:, 0]))
        bad = cr < 1e-12
        changed = bad.any()
        P = P[~bad] if changed else P
    return P


def area2(P):
    x, y = P[:, 0], P[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def earclip(P):
    """Triangulate a simple polygon (n, 2); returns index triples, or None.
    Triangles run along increasing index if area2(P) > 0, else decreasing."""
    n = len(P)
    idx = list(range(n)) if area2(P) > 0 else list(range(n))[::-1]
    tris, guard, i = [], 0, 0
    while len(idx) > 3:
        m = len(idx)
        if guard > m:
            return None
        a, b, c = idx[(i - 1) % m], idx[i % m], idx[(i + 1) % m]
        pa, pb, pc = P[a], P[b], P[c]
        cr = (pb[0] - pa[0]) * (pc[1] - pb[1]) - (pb[1] - pa[1]) * (pc[0] - pb[0])
        ear = cr > 1e-14
        if ear:
            others = np.array([j for j in idx if j not in (a, b, c)])
            if len(others):
                Q = P[others]
                d1 = (pb[0] - pa[0]) * (Q[:, 1] - pa[1]) - (pb[1] - pa[1]) * (Q[:, 0] - pa[0])
                d2 = (pc[0] - pb[0]) * (Q[:, 1] - pb[1]) - (pc[1] - pb[1]) * (Q[:, 0] - pb[0])
                d3 = (pa[0] - pc[0]) * (Q[:, 1] - pc[1]) - (pa[1] - pc[1]) * (Q[:, 0] - pc[0])
                ear = not np.any((d1 >= 0) & (d2 >= 0) & (d3 >= 0))
        if ear:
            tris.append((a, b, c))
            idx.pop(i % m)
            guard = 0
            i = i % max(len(idx), 1)
        else:
            i += 1
            guard += 1
    tris.append(tuple(idx))
    return tris


def section_loops(V, T, c, eu, ev, nrm):
    """Cross-section loops of a welded mesh in the plane (c, nrm), as 2D
    polygons in (eu, ev) coordinates about c. Drops sliver loops."""
    out, dropped = [], 0
    for L3 in slice_loops(V, T, c, nrm):
        Q = clean_loop(np.column_stack([(L3 - c) @ eu, (L3 - c) @ ev]))
        if len(Q) < 3 or abs(area2(Q)) < 1e-4:
            dropped += 1
            continue
        out.append(Q)
    return out, dropped


def loops_mask(loops, U, W):
    """Parity fill of (possibly nested) loops on grid points (U, W)."""
    from matplotlib.path import Path
    pts = np.column_stack([U.ravel(), W.ravel()])
    m = np.zeros(len(pts), bool)
    for q in loops:
        m ^= Path(q).contains_points(pts)
    return m.reshape(U.shape)


def loops_centroid(loops, h=0.05):
    allp = np.vstack(loops)
    lo, hi = allp.min(0) - h, allp.max(0) + h
    U, W = np.meshgrid(np.arange(lo[0], hi[0], h), np.arange(lo[1], hi[1], h), indexing="ij")
    m = loops_mask(loops, U, W)
    return np.array([U[m].mean(), W[m].mean()])


def loops_dice(A, B, h=0.05):
    if not A or not B:
        return 0.0
    allp = np.vstack(A + B)
    lo, hi = allp.min(0) - h, allp.max(0) + h
    U, W = np.meshgrid(np.arange(lo[0], hi[0], h), np.arange(lo[1], hi[1], h), indexing="ij")
    return dice(loops_mask(A, U, W), loops_mask(B, U, W))


def own_axis(V, T, z_top, span=30.0, step=5.0):
    """Straight axis of a tube near its top: centroids of horizontal sections
    from z_top - span to z_top - step, fitted with a line. Returns (point,
    unit direction pointing up)."""
    ex, ey, ez = np.eye(3)
    cs = []
    for z in np.arange(z_top - span, z_top - step + 1e-9, step):
        c0 = np.array([0.0, 0.0, z])
        loops, _ = section_loops(V, T, c0, ex, ey, ez)
        if not loops:
            continue
        allp = np.vstack(loops)
        h = 0.1
        U, W = np.meshgrid(np.arange(allp[:, 0].min(), allp[:, 0].max(), h),
                           np.arange(allp[:, 1].min(), allp[:, 1].max(), h), indexing="ij")
        m = loops_mask(loops, U, W)
        cs.append([U[m].mean(), W[m].mean(), z])
    cs = np.array(cs)
    p = cs.mean(0)
    d = np.linalg.svd(cs - p)[2][0]
    return p, d if d[2] > 0 else -d


def clip_indexed(V, T, p0, nrm):
    """The part of an indexed surface on the side (x - p0) . nrm < 0, crossing
    triangles cut exactly (RADO's swept STLs use triangles tens of mm long).
    Each cut point is computed once per mesh edge and shared by both faces of
    that edge. Returns (V', T', cut_edges): cut_edges are the new boundary
    edges (a, b) in the direction their own face runs them."""
    d = (V - p0) @ nrm
    d = np.where(np.abs(d) < 1e-9, 1e-9, d)
    below = d < 0
    sb = below[T]
    out = [T[sb.all(1)]]
    newV, cut, edges = [], {}, []

    def cp(u, w):
        key = (min(u, w), max(u, w))
        if key not in cut:
            a, b = key
            cut[key] = len(V) + len(newV)
            newV.append(V[a] + d[a] / (d[a] - d[b]) * (V[b] - V[a]))
        return cut[key]
    for tri in T[sb.any(1) & ~sb.all(1)]:
        poly, cuts = [], []
        for i in range(3):
            u, w = tri[i], tri[(i + 1) % 3]
            if below[u]:
                poly.append(u)
            if below[u] != below[w]:
                cuts.append(len(poly))
                poly.append(cp(u, w))
        out.append(np.array([[poly[0], poly[j], poly[j + 1]] for j in range(1, len(poly) - 1)]))
        i0, i1 = cuts
        # the boundary edge joins the two cut points, consecutive in poly order
        edges.append((poly[i0], poly[i1]) if (i1 - i0) % len(poly) == 1 else (poly[i1], poly[i0]))
    V2 = np.vstack([V, np.array(newV).reshape(-1, 3)])
    T2 = np.vstack(out)
    used = np.unique(T2)                       # drop vertices above the plane
    remap = np.full(len(V2), -1)
    remap[used] = np.arange(len(used))
    return V2[used], remap[T2], [(remap[a], remap[b]) for a, b in edges]


def chain_loops(edges):
    """Directed edges -> closed vertex loops in that direction. Raises if the
    cut is not a set of closed loops (an open source mesh)."""
    nxt = dict(edges)
    assert len(nxt) == len(edges), "cut edges branch: the source surface is not manifold here"
    loops, seen = [], set()
    for start in nxt:
        if start in seen:
            continue
        loop, cur = [], start
        while cur not in seen:
            seen.add(cur)
            loop.append(cur)
            assert cur in nxt, "cut loop is open: the source surface has a hole at the join"
            cur = nxt[cur]
        assert cur == start, "cut loops merge"
        loops.append(np.array(loop))
    return loops


def extend_body(tri, p0, eu, ev, nrm, step, t_end):
    """One closed surface: tri cut at the plane (p0, nrm), walls swept from the
    cut points along S^t for t in [0, t_end], and a cap at t_end. nrm points
    into the swept part. Returns (V, T, loops_uv, cap_methods)."""
    V, T = stlio.weld(tri)
    Vc, Tc, edges = clip_indexed(V, T, p0, nrm)
    loops = chain_loops(edges)
    ang = step.angle_deg() * t_end
    ts = np.linspace(0, t_end, max(2, int(np.ceil(ang / SWEEP_MAX_DEG)) + 1))
    Ms = [step.power(t) for t in ts[1:]]
    verts, faces, nv = [Vc], [Tc], len(Vc)
    loops_uv, methods = [], []
    for L in loops:
        ring0 = Vc[L]
        n = len(L)
        rings = [L] + [nv + i * n + np.arange(n) for i in range(len(Ms))]
        verts.append(np.vstack([apply(M, ring0) for M in Ms]))
        nv += n * len(Ms)
        j, k = np.arange(n), (np.arange(n) + 1) % n
        for r0, r1 in zip(rings[:-1], rings[1:]):
            # The cut face runs j -> k, so the wall must run k -> j.
            faces.append(np.column_stack([r0[k], r0[j], r1[j]]))
            faces.append(np.column_stack([r0[k], r1[j], r1[k]]))
        Q = np.column_stack([(ring0 - p0) @ eu, (ring0 - p0) @ ev])
        loops_uv.append(Q)
        top = rings[-1]
        ct = earclip(Q)
        if ct is not None:
            ct, flip = np.array(ct), area2(Q) > 0
            methods.append("earclip")
        else:
            # Centroid fan: still right for the ray-PARITY containment test
            # fem/ uses (each interior point is covered an odd number of times).
            verts.append(apply(Ms[-1], ring0.mean(0)[None]))
            top = np.concatenate([top, [nv]])
            nv += 1
            ct, flip = np.array([(n, a, (a + 1) % n) for a in range(n)]), True
            methods.append("fan")
        # The walls end running j -> k on the top ring, so the cap runs k -> j.
        faces.append(top[ct[:, ::-1]] if flip else top[ct])
    return np.vstack(verts), np.vstack(faces), loops_uv, methods


def seg_axis(P):
    """(centroid, unit axis pointing up, axial coordinates) of a thin tube."""
    c = P.mean(0)
    d = np.linalg.svd(P - c, full_matrices=False)[2][0]
    d = d if d[2] > 0 else -d
    return c, d, (P - c) @ d


def seg_ends(P, band=0.3):
    c, d, s = seg_axis(P)
    return P[s < s.min() + band].mean(0), P[s > s.max() - band].mean(0), d


def join_trunk(tri, ganglion_centre, step):
    """Bend the upper end of a sympathetic trunk segment so that it overlaps
    its own S-copy (the segment of the next level up) by CHAIN_OVERLAP_MM.
    The displacement is w(u) * g, w a smoothstep in the axial coordinate
    from the ganglion (0) to the top end (1): the part below the ganglion,
    the ganglion's own stretch and the bottom end do not move."""
    P = tri.reshape(-1, 3)
    c, d, s = seg_axis(P)
    _, top, _ = seg_ends(P)
    bottom_up, _, d_up = seg_ends(apply(step.M, P))
    g = bottom_up + CHAIN_OVERLAP_MM * d_up - top
    s_g = (ganglion_centre - c) @ d
    u = np.clip((s - s_g) / (s.max() - s_g), 0, 1)
    w = u * u * (3 - 2 * u)
    Pj = P + w[:, None] * g
    bend = float(np.degrees(np.arctan(1.5 * np.linalg.norm(g - (g @ d) * d) / (s.max() - s_g))))
    return Pj.reshape(-1, 3, 3), {"gap_closed_mm": round(float(np.linalg.norm(g)), 3),
                                  "axial_mm": round(float(g @ d), 3),
                                  "lateral_mm": round(float(np.linalg.norm(g - (g @ d) * d)), 3),
                                  "bent_length_mm": round(float(s.max() - s_g), 2),
                                  "max_bend_deg": round(bend, 2)}


# --------------------------------------------------------------------------
# Checks
# --------------------------------------------------------------------------
class TesterCache:
    def __init__(self):
        self._c = {}

    def get(self, path):
        if path not in self._c:
            self._c[path] = InsideTester(*stlio.read_stl(path))
        return self._c[path]


def section_grid(c, eu, ev, half_u, v_lo, v_hi, h):
    u = np.arange(-half_u, half_u + h / 2, h)
    v = np.arange(v_lo, v_hi + h / 2, h)
    U, W = np.meshgrid(u, v, indexing="ij")
    return c + U.ravel()[:, None] * eu + W.ravel()[:, None] * ev, U.shape


def label_points(pts, bodies_by_pipe, cache, order):
    lab = np.full(len(pts), -1, dtype=np.int16)
    todo = np.arange(len(pts))
    lo, hi = pts.min(0), pts.max(0)
    for li, pipe in enumerate(order):
        if not todo.size:
            break
        m = np.zeros(todo.size, bool)
        for path, bb in bodies_by_pipe.get(pipe, []):
            if np.any(bb[1] < lo) or np.any(bb[0] > hi):
                continue
            m |= cache.get(path)(pts[todo])
        lab[todo[m]] = li
        todo = todo[~m]
    return lab


def dice(a, b):
    s = a.sum() + b.sum()
    return float(2 * (a & b).sum() / s) if s else float("nan")


def disc_contact(disc_path, below_path, above_path, cache, probe=1.5, dprobe=0.05):
    """How a disc sits between its vertebrae. For each face of the disc
    facing a vertebra (|normal . disc axis| > 0.7), the distance along its
    outward normal to that vertebra's bone: 0 = touching or embedded.
    Returned as area-weighted fractions and percentiles."""
    V, T = stlio.read_stl(disc_path)
    a, b, c = V[T[:, 0]], V[T[:, 1]], V[T[:, 2]]
    n = np.cross(b - a, c - a)
    area = 0.5 * np.linalg.norm(n, axis=1)
    n = n / np.maximum(2 * area, 1e-30)[:, None]
    cen = (a + b + c) / 3
    axis = np.linalg.svd(V - V.mean(0), full_matrices=False)[2][2]
    if (V[np.argmax(V[:, 2])] - V.mean(0)) @ axis < 0:
        axis = -axis
    out = {}
    for side, path, sgn in (("lower", below_path, -1), ("upper", above_path, 1)):
        face = (n @ axis) * sgn > 0.7
        if path is None or not face.any():
            out[side] = None
            continue
        tester = cache.get(path)
        gap = np.full(face.sum(), np.inf)
        for dd in np.arange(0, probe + 1e-9, dprobe):
            hit = tester(cen[face] + dd * n[face])
            gap = np.where(np.isinf(gap) & hit, dd, gap)
        w = area[face] / area[face].sum()
        srt = np.argsort(gap)
        cw = np.cumsum(w[srt])
        out[side] = {"vertebra": os.path.basename(path),
                     "facing_area_mm2": round(float(area[face].sum()), 1),
                     "touching_fraction": round(float(w[gap <= dprobe].sum()), 4),
                     "gap_p50_mm": float(gap[srt][np.searchsorted(cw, 0.50)]),
                     "gap_p95_mm": float(gap[srt][np.searchsorted(cw, 0.95)]),
                     "beyond_probe_fraction": round(float(w[np.isinf(gap)].sum()), 4)}
    return out


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--levels", type=int, default=4, help="vertebral levels to add above the top")
    ap.add_argument("--join-back", type=float, default=10.0,
                    help="mm below the lowest top of the continuous bodies where all of them "
                         "are cut and the sweeps start. RADO's continuous bodies are not the "
                         "constant section within ~3-8 mm of their ends (measured), so keep "
                         "this >= 10")
    ap.add_argument("--out", default=None)
    ap.add_argument("--no-check", action="store_true", help="skip the section checks and figure")
    args = ap.parse_args()
    t0 = time.time()
    K = args.levels
    warnings = []
    out = args.out or os.path.join(C.OUT, "geometry", "rado_up%d" % K)
    if os.path.exists(out):
        shutil.rmtree(out)
    os.makedirs(out)
    entries = parse_tissue_map(C.TISSUE_MAP_YAML)

    # ---- load every source STL -------------------------------------------
    bodies = []
    for fn in sorted(os.listdir(C.STL)):
        if not fn.lower().endswith(".stl"):
            continue
        path = os.path.join(C.STL, fn)
        stem = os.path.splitext(fn)[0]
        tri = read_raw(path)
        cls = class_of(stem, entries)
        bodies.append({"file": fn, "path": path, "stem": stem, "cls": cls,
                       "pipe": PIPE_OF_CLASS.get(cls), "tri": tri,
                       "centroid": tri.reshape(-1, 3).mean(0),
                       "zspan": float(np.ptp(tri[:, :, 2]))})
    anatomy = [b for b in bodies if b["pipe"]]
    continuous = [b for b in anatomy if b["zspan"] > CONTINUOUS_SPAN_MM]
    per_level = [b for b in anatomy if b["zspan"] <= CONTINUOUS_SPAN_MM]
    print("%d source STLs: %d anatomical (%d continuous, %d per-level), %d lead hardware/other"
          % (len(bodies), len(anatomy), len(continuous), len(per_level), len(bodies) - len(anatomy)))

    # ---- the vertebral step S ------------------------------------------------
    verts = sorted([b for b in per_level if b["pipe"] == "vertebra"], key=lambda b: b["centroid"][2])
    P = [b["tri"].reshape(-1, 3) for b in verts]
    rado_steps = []
    for a, b, va, vb in zip(P[:-1], P[1:], verts[:-1], verts[1:]):
        R, t, res = rigid_fit(a, b)
        assert res < RIGID_TOL_MM, "vertebrae are no longer exact rigid copies"
        M = np.eye(4)
        M[:3, :3], M[:3, 3] = R, t
        st = Step(M)
        rado_steps.append({"from": va["stem"], "to": vb["stem"], "matrix": M.round(9).tolist(),
                           "angle_deg": round(st.angle_deg(), 4), "axis": st.axis().round(5).tolist(),
                           "centre_yz_mm": st.centre_yz().round(2).tolist(),
                           "radius_at_vertebra_centroid_mm": round(float(np.linalg.norm(
                               vb["centroid"][1:] - st.centre_yz())), 1),
                           "fit_residual_max_mm": res})
        print("  RADO step %s -> %s: %.3f deg about %s, centre (y, z) = (%.1f, %.1f), radius %.0f mm at the vertebra"
              % (va["stem"][8:], vb["stem"][8:], rado_steps[-1]["angle_deg"], rado_steps[-1]["axis"],
                 *rado_steps[-1]["centre_yz_mm"], rado_steps[-1]["radius_at_vertebra_centroid_mm"]))
    step = st          # the top step, T11 -> T10
    print("new levels are S^k, S = the top step")

    # ---- centreline (selection, checks and figure only) ----------------------
    cpaths = C.TISSUE_BODIES["white"] + C.TISSUE_BODIES["grey"]
    cen = cord_centroids(cpaths)
    coef = np.polyfit(cen[:, 2], cen[:, 1], 3)
    fit_res = cen[:, 1] - np.polyval(coef, cen[:, 2])
    x0 = float(cen[:, 0].mean())
    # The cubic runs past the cord data to the top of the tallest continuous
    # body (the aorta, 12 mm above the cord), so every body's top projects to
    # its own arc length instead of clamping at the end of the data.
    z_hi = max(float(b["tri"][:, :, 2].max()) for b in continuous) + 5.0
    curve = Curve.cubic(coef, x0, cen[0, 2] - 12.0, z_hi)
    print("centreline: cubic y(z) residual rms %.4f mm max %.4f mm; x drift sd %.4f mm"
          % (np.sqrt(np.mean(fit_res ** 2)), np.abs(fit_res).max(), cen[:, 0].std()))
    s_vref = float(curve.project(verts[-1]["centroid"])[0])
    D = float(s_vref - curve.project(verts[-2]["centroid"])[0])

    # ---- group per-level bodies into families of rigid copies -------------
    fams = []          # list of lists of bodies
    for b in sorted(per_level, key=lambda b: b["stem"]):
        if b["pipe"] == "disc":
            continue
        A = b["tri"].reshape(-1, 3)
        for f in fams:
            r = f[0]["tri"].reshape(-1, 3)
            if r.shape == A.shape and rigid_fit(r, A)[2] < RIGID_TOL_MM:
                f.append(b)
                break
        else:
            fams.append([b])
    groups = {}
    for f in fams:
        if len(f) > 1:
            groups["family:" + f[0]["stem"]] = f
        else:
            key = "name:" + re.sub(r"-\d+(?= |$)", "-#", f[0]["stem"], count=1)
            groups.setdefault(key, []).extend(f)
    groups["discs"] = [b for b in per_level if b["pipe"] == "disc"]

    # Which members form the top level: measured from the top vertebra's
    # centroid along the centreline. RADO's levels are rigid copies but were
    # NOT placed by one repeated transform (e.g. the left DRG steps 31.3 /
    # 27.6 / 25.8 mm with sideways shifts while the vertebrae step 31.9 /
    # 30.8 mm), so the TOP level is copied as one rigid snapshot of RADO's own
    # top-level arrangement.
    copy_set, group_report = [], []
    for key, mem in sorted(groups.items()):
        sc = np.array([curve.project(b["centroid"])[0] for b in mem])
        top = sc.max()
        offset = (top - s_vref) / D
        tops = [b for b, s in zip(mem, sc) if s > top - 0.5 * D]
        include = offset > -0.5
        if -0.65 < offset < -0.35:
            warnings.append("group %s top member is ambiguously placed (%.2f levels from the top vertebra)" % (key, offset))
        group_report.append({"group": key, "class": mem[0]["pipe"], "members": len(mem),
                             "top_members": [b["stem"] for b in tops],
                             "top_offset_levels_from_top_vertebra_centroid": round(float(offset), 3),
                             "copied": bool(include)})
        if include:
            copy_set.extend(tops)
    print("per-level groups: %d (%d copied as the top-level unit, %d bodies per new level)"
          % (len(groups), sum(g["copied"] for g in group_report), len(copy_set)))

    # ---- sympathetic trunk joins ------------------------------------------
    trunks = [b for b in copy_set if "sympathetic_a" in b["stem"]]
    ganglia = [b for b in copy_set if "sympathetic_ball" in b["stem"]]
    joined, chain_report = {}, []
    for b in trunks:
        gang = min(ganglia, key=lambda g: np.linalg.norm(g["centroid"] - b["centroid"]))
        tri_j, info = join_trunk(b["tri"], gang["centroid"], step)
        joined[b["stem"]] = tri_j
        info.update({"segment": b["stem"], "ganglion": gang["stem"]})
        chain_report.append(info)
        print("  trunk %s: top end moved %.2f mm (axial %.2f, lateral %.2f) over %.1f mm above its ganglion, "
              "bend <= %.1f deg" % (b["stem"][8:], info["gap_closed_mm"], info["axial_mm"], info["lateral_mm"],
                                    info["bent_length_mm"], info["max_bend_deg"]))

    # ---- level names -------------------------------------------------------
    top_disc = max(groups["discs"], key=lambda b: b["centroid"][2])
    m = re.search(r"T_disk_(\d+)_(\d+|L)", top_disc["stem"])
    top_vert_disc_naming = int(m.group(2))          # the vertebra below the top disc
    levels = []
    for k in range(1, K + 1):
        v = top_vert_disc_naming - k
        Mk = step.power(k)
        levels.append({"k": k, "vertebra_disc_naming": "T%d" % v, "vertebra_paper_naming": "T%d" % (v - 1),
                       "disc_above_disc_naming": "T%d/T%d" % (v - 1, v),
                       "vertebra_centroid_mm": apply(Mk, verts[-1]["centroid"][None])[0].round(3).tolist(),
                       "transform_from_top_level": Mk.round(9).tolist(),
                       "rot_deg_from_top_level": round(k * step.angle_deg(), 3)})

    # ---- write: RADO's files, minus the ones replaced ----------------------
    replaced = {b["file"] for b in continuous} | {b["file"] for b in trunks}
    outputs = []
    for b in bodies:
        if b["file"] not in replaced:
            shutil.copy2(b["path"], os.path.join(out, b["file"]))
    for b in trunks:
        stem = b["stem"] + "__joined"
        assert class_of(stem, entries) == b["cls"], stem
        write_stl(os.path.join(out, stem + ".STL"), joined[b["stem"]])
        outputs.append({"file": stem + ".STL", "kind": "joined", "source": b["file"], "class": b["pipe"],
                        "triangles": len(joined[b["stem"]])})
    for L in levels:
        k = L["k"]
        M = np.array(L["transform_from_top_level"])
        for b in copy_set:
            if b["pipe"] == "disc":
                u, lo = re.search(r"T_disk_(\d+)_(\d+|L)", b["stem"]).groups()
                lo = int(lo) if lo != "L" else 13
                stem = "T8-10 - T_disk_%d_%dx-1__up%d" % (int(u) - k, lo - k, k)
            else:
                stem = "%s__up%d" % (b["stem"], k)
            assert class_of(stem, entries) == b["cls"], stem
            # every trunk copy but the topmost carries the bend that joins it upward
            src = joined[b["stem"]] if b["stem"] in joined and k < K else b["tri"]
            tri = apply(M, src.reshape(-1, 3)).reshape(-1, 3, 3)
            write_stl(os.path.join(out, stem + ".STL"), tri)
            outputs.append({"file": stem + ".STL", "kind": "copy", "level_k": k, "source": b["file"],
                            "class": b["pipe"], "triangles": len(tri)})
    print("wrote %d level copies" % sum(o["kind"] == "copy" for o in outputs))

    # ---- extend the continuous bodies ---------------------------------------
    # Canal contents are cut perpendicular to the cord centreline; a body that
    # is not constant in that frame (the aorta runs on its own path) is cut
    # perpendicular to its own top axis instead. Measured per body. Either way
    # the sweep follows the S-orbit, as the vertebrae do.
    # All are cut at ONE plane, --join-back below the lowest top: the nested
    # canal layers (cord, CSF, dura, epidural) then keep sections taken at the
    # same place, which are consistent with each other. Cut at their own tops
    # (up to 1.8 mm apart) they left 2.3x more unassigned canal points. Each
    # body still ends where RADO's own does relative to the top level (the
    # aorta 12 mm above the cord).
    s_tops = {b["file"]: float(curve.project(b["tri"].reshape(-1, 3)).max()) for b in continuous}
    s_cut = min(s_tops.values()) - args.join_back
    ext_report = []
    for b in continuous:
        V, T = stlio.weld(b["tri"])
        s_top = s_tops[b["file"]]
        s_join = s_cut
        ex = np.array([1.0, 0, 0])
        loops_hi, _ = section_loops(V, T, curve.point(s_join), ex, curve.normal(s_join), curve.tangent(s_join))
        loops_lo, _ = section_loops(V, T, curve.point(s_join - 20), ex, curve.normal(s_join - 20),
                                    curve.tangent(s_join - 20))
        constancy = loops_dice(loops_hi, loops_lo)
        drift = float(np.linalg.norm(loops_centroid(loops_hi) - loops_centroid(loops_lo)))
        if drift < 1.0:
            method = "centreline"
            p0, eu, ev, nrm = curve.point(s_join), ex, curve.normal(s_join), curve.tangent(s_join)
        else:
            method = "own_axis"
            z_top = float(V[:, 2].max())
            p_axis, nrm = own_axis(V, T, z_top)
            p0 = p_axis + nrm * ((z_top - args.join_back - p_axis[2]) / nrm[2])
            eu = ex - (ex @ nrm) * nrm
            eu /= np.linalg.norm(eu)
            ev = np.cross(nrm, eu)
        orbit = step.L[:3, :3] @ p0 + step.L[:3, 3]
        turn = float(np.degrees(np.arccos(np.clip(orbit @ nrm / np.linalg.norm(orbit), -1, 1))))
        below_top = s_top - s_join if method == "centreline" else args.join_back
        t_end = K + below_top / step.speed(p0)
        Vx, Tx, loops_uv, methods = extend_body(b["tri"], p0, eu, ev, nrm, step, t_end)
        stem = "%s__extended_up%d" % (b["stem"], K)
        assert class_of(stem, entries) == b["cls"], stem
        write_stl(os.path.join(out, stem + ".STL"), Vx[Tx])
        topo = topology(Vx, Tx)
        topo_file = topology(*stlio.weld(read_raw(os.path.join(out, stem + ".STL")), tol=1e-5))
        topo_src = topology(V, T)
        ext_report.append({"source": b["file"], "class": b["pipe"], "file": stem + ".STL", "method": method,
                           "constancy_dice_in_centreline_frame_20mm": round(constancy, 4),
                           "section_centroid_drift_in_centreline_frame_20mm": round(drift, 4),
                           "turn_at_join_deg": round(turn, 3),
                           "join_frame": {"centre": p0.tolist(), "e_u": eu.tolist(), "e_v": ev.tolist(),
                                          "normal": nrm.tolist()},
                           "cut_below_own_top_mm": round(below_top, 3), "t_end_levels": round(t_end, 4), "mm_per_level_at_join": round(step.speed(p0), 3),
                           "loops": len(loops_uv), "loop_points": [len(q) for q in loops_uv], "caps": methods,
                           "topology": topo, "topology_as_written_float32": topo_file,
                           "topology_source": topo_src,
                           "loops_uv": [q.round(5).tolist() for q in loops_uv] if b["pipe"] == "epidural" else None})
        outputs.append({"file": stem + ".STL", "kind": "extended", "source": b["file"], "class": b["pipe"],
                        "triangles": len(Tx)})
        bad = {k: topo_file[k] - topo_src[k] for k in ("boundary_edges", "nonmanifold_edges", "misoriented_edges")}
        if any(bad.values()):
            warnings.append("%s: not a clean closed surface as written (%s more than the source)" % (stem, bad))
        if "fan" in methods:
            warnings.append("%s: a cap used the fan fallback (parity-correct, geometrically folded)" % stem)
        print("  %-10s %-28s cut %.1f mm below its top, turns %.2f deg there; %d loop(s); closed: %d open / "
              "%d non-manifold / %d misoriented edges (source %d/%d/%d)"
              % (method, b["stem"][8:], below_top, turn, len(loops_uv), topo_file["boundary_edges"],
                 topo_file["nonmanifold_edges"], topo_file["misoriented_edges"], topo_src["boundary_edges"],
                 topo_src["nonmanifold_edges"], topo_src["misoriented_edges"]))

    for o in outputs:
        o["sha256"] = sha256(os.path.join(out, o["file"]))

    # The checks' and figure's centreline: the cubic, then the S-orbit of its
    # point at the cord's join.
    wm = next(r for r in ext_report if r["class"] == "white")
    s_join_wm = float(curve.project(np.array(wm["join_frame"]["centre"]))[0])
    curve_ext = curve.continued(s_join_wm, step, wm["t_end_levels"])
    write_centreline(curve_ext, s_join_wm, os.path.join(out, "centreline.json"))

    # ---- checks ------------------------------------------------------------
    checks = {}
    if not args.no_check:
        checks = run_checks(curve_ext, step, out, bodies, replaced, outputs, levels, s_vref, D,
                            ext_report, chain_report, warnings)

    manifest = {
        "tool": "src/extend/extend_rado.py",
        "created": time.strftime("%Y-%m-%d %H:%M:%S"),
        "parameters": {"levels": K, "join_back_mm": args.join_back, "chain_overlap_mm": CHAIN_OVERLAP_MM,
                       "sweep_max_deg_between_rings": SWEEP_MAX_DEG},
        "inputs": {b["file"]: sha256(b["path"]) for b in bodies},
        "tissue_map_sha256": sha256(C.TISSUE_MAP_YAML),
        "step": {"note": "new level k = S^k applied to RADO's top-level unit; S = RADO's top vertebral step",
                 "rado_steps": rado_steps, "S": step.M.round(9).tolist(), "S_angle_deg": step.angle_deg()},
        "centreline": {"model": "x = x0; y = cubic(z) inside the cord, then the S-orbit (checks/figure only)",
                       "x0_mm": x0, "cubic_coef_high_first": coef.tolist(),
                       "fit_residual_rms_mm": float(np.sqrt(np.mean(fit_res ** 2))),
                       "fit_residual_max_mm": float(np.abs(fit_res).max()),
                       "x_drift_sd_mm": float(cen[:, 0].std()),
                       "data_z_range_mm": [float(cen[0, 2]), float(cen[-1, 2])]},
        "level_spacing_arc_mm": D,
        "new_levels": levels,
        "level_naming_note": "disc naming (default, used in new disc filenames) makes RADO's "
                             "vertebrae T10-T12; the paper's naming is one level higher. "
                             "See src/freecad/README.md 'Vertebral levels'.",
        "per_level_groups": group_report,
        "replaced_rado_files": sorted(replaced),
        "sympathetic_joins": chain_report,
        "continuous_extensions": ext_report,
        "outputs": outputs,
        "checks": checks,
        "warnings": warnings,
        "elapsed_s": round(time.time() - t0, 1),
    }
    with open(os.path.join(out, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1)
    print("\n%d warning(s)" % len(warnings))
    for w in warnings:
        print("  - " + w)
    print("wrote %s (%d new STLs, %d RADO files replaced) in %.0f s" % (out, len(outputs), len(replaced),
                                                                       time.time() - t0))


def write_centreline(curve, s_join, path, ds=0.5):
    """The cord centreline of the extended set, sampled for other tools: the
    cubic fit inside RADO's cord, then the S-orbit of its point at the cut."""
    s = np.arange(curve.S[0], curve.S[-1], ds)
    P = np.array([curve.point(x) for x in s])
    T = np.array([curve.tangent(x) for x in s])
    with open(path, "w") as f:
        json.dump({"frame": "RADO STL mm, LPS: +x patient left, +y posterior, +z superior",
                   "definition": "x = x0; y = cubic(z) up to arc length s_cut, then the orbit of that point "
                                 "under S^t (manifest.json step.S); tangent points superior",
                   "x0_mm": curve.x0, "s_cut_mm": round(float(s_join), 4), "ds_mm": ds,
                   "s_mm": s.round(4).tolist(), "points_mm": P.round(5).tolist(),
                   "tangents": T.round(6).tolist()}, f)


def run_checks(curve, step, out, bodies, replaced, outputs, levels, s_vref, D, ext_report, chain_report, warnings):
    """Each check aimed at one way the extension could be wrong.

    1. joins: each extended body against the original, just above the cut
       (where the original still exists) and just below it (must be identical).
    2. copies: RADO's top unit + the extended bodies, sectioned and carried by
       S^k, against copy k + the extended bodies on the same points.
    3. sympathetic trunk: every junction overlaps.
    4. discs: how each disc meets its vertebrae, RADO's own against the new.
    5. gaps: in the new region, every point inside the canal's outer outline
       must belong to some tissue (a gap would become 0.04 S/m background).
    """
    entries = parse_tissue_map(C.TISSUE_MAP_YAML)
    order = [p for p, _ in C.ANATOMY_CLASS]
    cache = TesterCache()
    files = [b["file"] for b in bodies if b["file"] not in replaced] + [o["file"] for o in outputs]

    def index(fs):
        by = {}
        for fn in fs:
            path = os.path.join(out, fn)
            pipe = PIPE_OF_CLASS.get(class_of(os.path.splitext(fn)[0], entries))
            if not pipe:
                continue
            v = stlio.read_stl(path)[0]
            by.setdefault(pipe, []).append((path, (v.min(0), v.max(0))))
        return by
    res = {}

    # 1. joins -----------------------------------------------------------------
    joins = []
    for r in ext_report:
        f = r["join_frame"]
        c, eu, ev, nrm = (np.array(f[k]) for k in ("centre", "e_u", "e_v", "normal"))
        src, ext = os.path.join(C.STL, r["source"]), os.path.join(out, r["file"])
        row = {"body": r["source"]}
        for name, dz in (("below_cut", -5.0), ("above_cut", 5.0)):
            cc = c + dz * nrm
            V, T = stlio.read_stl(src)
            loops, _ = section_loops(V, T, cc, eu, ev, nrm)
            allp = np.vstack(loops)
            h = 0.05
            u = np.arange(allp[:, 0].min() - .5, allp[:, 0].max() + .5, h)
            w = np.arange(allp[:, 1].min() - .5, allp[:, 1].max() + .5, h)
            U, W = np.meshgrid(u, w, indexing="ij")
            pts = cc + U.ravel()[:, None] * eu + W.ravel()[:, None] * ev
            a, b = cache.get(src)(pts), cache.get(ext)(pts)
            row[name] = {"dice": round(dice(a, b), 4), "xor_area_mm2": round(float((a ^ b).sum() * h * h), 4)}
        joins.append(row)
    res["joins"] = joins
    print("joins, extended vs original 5 mm below / above the cut (Dice): " + ", ".join(
        "%s %.3f/%.3f" % (j["body"][8:].split("-")[0][:13], j["below_cut"]["dice"], j["above_cut"]["dice"])
        for j in joins))
    for j in joins:
        if j["below_cut"]["dice"] < 0.999:
            warnings.append("%s: extended body differs from the original below the cut" % j["body"])

    # 2. copies ------------------------------------------------------------------
    ext_files = [r["file"] for r in ext_report]
    unit_src = sorted({o["source"] for o in outputs if o["kind"] == "copy"})
    unit_src = [s if s not in replaced else s[:-4] + "__joined.STL" for s in unit_src]
    by_ref = index(unit_src + ext_files)
    per_level_cls = [order.index(c) for c in ("drg", "root", "sympathetic", "disc", "vertebra")]
    # The extended bodies end --join-back mm above the topmost unit, so the
    # topmost level's sections beyond that compare against an empty end of
    # the model (and its trunk copy is deliberately left unbent). Recorded,
    # not gated.
    f_top = min((r["t_end_levels"] - len(levels)) * r["mm_per_level_at_join"] for r in ext_report) / D
    reps = []
    for L in levels:
        M = np.array(L["transform_from_top_level"])
        by_k = index([o["file"] for o in outputs if o["kind"] == "copy" and o["level_k"] == L["k"]] + ext_files)
        for f in (-0.25, 0.0, 0.25, 0.5, 0.75):
            s = s_vref + f * D
            p0, _ = section_grid(curve.point(s), np.array([1.0, 0, 0]), curve.normal(s), 45, -45, 40, 0.2)
            l0 = label_points(p0, by_ref, cache, order)
            l1 = label_points(apply(M, p0), by_k, cache, order)
            occ = (l0 >= 0) | (l1 >= 0)
            pl = np.isin(l0, per_level_cls)
            reps.append({"level_k": L["k"], "offset_levels": f,
                         "agreement_where_occupied": round(float((l0 == l1)[occ].mean()), 5),
                         "agreement_on_per_level_tissue": round(float((l0 == l1)[pl].mean()), 5)
                         if pl.any() else None, "points_occupied": int(occ.sum()),
                         "beyond_model_top": bool(L["k"] == len(levels) and f > f_top)})
        mine = [r for r in reps if r["level_k"] == L["k"] and not r["beyond_model_top"]]
        worst_pl = min(r["agreement_on_per_level_tissue"] for r in mine if r["agreement_on_per_level_tissue"] is not None)
        worst = min(r["agreement_where_occupied"] for r in mine)
        print("  copy check level +%d (%s): per-level tissue agreement >= %.4f; all tissue >= %.4f"
              % (L["k"], L["vertebra_disc_naming"], worst_pl, worst))
        if worst_pl < 0.98:
            warnings.append("copy check level +%d: per-level tissue agreement %.4f" % (L["k"], worst_pl))
    res["copies"] = reps

    # 3. sympathetic trunk ---------------------------------------------------------
    # Each junction: the lower segment's top end centre, pulled back 0.2 mm
    # along its axis, must lie inside the segment above.
    seg_files = {}
    for o in outputs:
        if "sympathetic_a" in o["file"]:
            src = o["source"]
            seg_files.setdefault(src, {})[o.get("level_k", 0)] = o["file"]
    junctions = []
    for src, by_k in sorted(seg_files.items()):
        for k in range(0, len(levels)):
            lo, hi = os.path.join(out, by_k[k]), os.path.join(out, by_k[k + 1])
            Pl = read_raw(lo).reshape(-1, 3)
            _, top, d = seg_ends(Pl)
            inside = bool(cache.get(hi)((top - 0.2 * d)[None])[0])
            dist = float(cKDTree(stlio.read_stl(hi)[0]).query(top)[0])
            junctions.append({"lower": by_k[k], "upper": by_k[k + 1], "top_end_inside_upper": inside})
            if not inside:
                warnings.append("sympathetic junction %s -> %s does not overlap (end %.2f mm away)"
                                % (by_k[k], by_k[k + 1], dist))
    res["sympathetic_junctions"] = junctions
    print("  sympathetic trunk: %d of %d junctions overlap" % (sum(j["top_end_inside_upper"] for j in junctions),
                                                              len(junctions)))

    # 4. discs ---------------------------------------------------------------------
    disc_files = [fn for fn in files if PIPE_OF_CLASS.get(class_of(os.path.splitext(fn)[0], entries)) == "disc"]
    vert_files = [fn for fn in files if PIPE_OF_CLASS.get(class_of(os.path.splitext(fn)[0], entries)) == "vertebra"]
    vc = {fn: read_raw(os.path.join(out, fn)).reshape(-1, 3).mean(0) for fn in vert_files}
    discs = []
    for fn in disc_files:
        dc = read_raw(os.path.join(out, fn)).reshape(-1, 3).mean(0)
        s_d = curve.project(dc)[0]
        below = [v for v in vert_files if curve.project(vc[v])[0] < s_d]
        above = [v for v in vert_files if curve.project(vc[v])[0] > s_d]
        b = max(below, key=lambda v: curve.project(vc[v])[0]) if below else None
        a = min(above, key=lambda v: curve.project(vc[v])[0]) if above else None
        r = disc_contact(os.path.join(out, fn), b and os.path.join(out, b), a and os.path.join(out, a), cache)
        r["disc"] = fn
        r["new"] = "__up" in fn
        discs.append(r)
    res["discs"] = discs
    print("  discs (fraction of facing area touching bone, lower/upper; 95th percentile gap mm):")
    for r in sorted(discs, key=lambda r: curve.project(read_raw(os.path.join(out, r["disc"])).reshape(-1, 3).mean(0))[0]):
        f = lambda s: "%.2f (p95 %.2f)" % (r[s]["touching_fraction"], r[s]["gap_p95_mm"]) if r[s] else "   -   "
        print("     %-34s %s  lower %s  upper %s" % (r["disc"][8:42], "new " if r["new"] else "RADO", f("lower"), f("upper")))

    # 5. gaps ----------------------------------------------------------------------
    # Sections along the epidural extension's own frames above its cut, and
    # along the cubic centreline below it; outline = the epidural outer loop.
    epi = next(r for r in ext_report if r["class"] == "epidural")
    loops = [np.array(q) for q in epi["loops_uv"]]
    outer = max(loops, key=lambda q: abs(area2(q)))
    jf = epi["join_frame"]
    c0, eu0, ev0 = (np.array(jf[k]) for k in ("centre", "e_u", "e_v"))
    lo, hi = outer.min(0), outer.max(0)
    U, W = np.meshgrid(np.arange(lo[0], hi[0], 0.1), np.arange(lo[1], hi[1], 0.1), indexing="ij")
    inside = loops_mask([outer], U, W).ravel()
    uv = np.column_stack([U.ravel()[inside], W.ravel()[inside]])
    by_ext = index(files)
    gaps = []
    s_join = float(curve.project(c0)[0])
    frames = [("rado", curve.point(s), np.array([1.0, 0, 0]), curve.normal(s))
              for s in np.arange(s_vref - 2.0 * D, s_join - 2.0, 5.0)]
    for t in np.arange(0.1, epi["t_end_levels"] - 0.15, 5.0 / epi["mm_per_level_at_join"]):
        M = step.power(t)
        frames.append(("extension", apply(M, c0[None])[0], M[:3, :3] @ eu0, M[:3, :3] @ ev0))
    for region, c, eu, ev in frames:
        pts = c + uv[:, :1] * eu + uv[:, 1:] * ev
        lab = label_points(pts, by_ext, cache, order)
        gaps.append({"z_mm": round(float(c[2]), 2), "region": region, "canal_points": int(len(pts)),
                     "background_points": int((lab == -1).sum())})
    # Background inside the canal outline is NOT zero in RADO itself: its
    # compartments are not tessellated compatibly, so thin unassigned slivers
    # sit along shared interfaces and around the root sleeves. The test is
    # that the extension adds none beyond RADO's own baseline.
    rb = [g["background_points"] for g in gaps if g["region"] == "rado"]
    eb = [g["background_points"] for g in gaps if g["region"] == "extension"]
    res["gaps_summary"] = {"grid_mm": 0.1, "rado_sections": len(rb), "rado_mean": float(np.mean(rb)),
                           "rado_max": int(max(rb)), "extension_sections": len(eb),
                           "extension_mean": float(np.mean(eb)), "extension_max": int(max(eb))}
    print("  gap check (background points inside the canal outline per section, 0.1 mm grid): "
          "RADO's own canal mean %.1f max %d; extension mean %.1f max %d"
          % (np.mean(rb), max(rb), np.mean(eb), max(eb)))
    if np.mean(eb) > 1.5 * np.mean(rb) + 2 or max(eb) > max(rb) + 5:
        warnings.append("gap check: extension has more unassigned canal points than RADO's own baseline")
    res["gaps"] = gaps

    try:
        plot(curve, step, by_ext, cache, order, entries, levels, s_vref, out)
        res["figure"] = "extension_check.png"
    except Exception as e:      # the figure is a convenience, never a gate
        warnings.append("figure failed: %r" % (e,))
    return res


def plot(curve, step, by_ext, cache, order, entries, levels, s_vref, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    rgb = {}
    for t in entries:
        pipe = PIPE_OF_CLASS.get(t["name"])
        if pipe and t["color_rgb"]:
            rgb[pipe] = np.array(t["color_rgb"]) / 255.0
    cols = [(0.08, 0.08, 0.1)] + [tuple(rgb.get(p, (0.5, 0.5, 0.5))) for p in order]
    cmap = ListedColormap(cols)
    h = 0.4
    # left: the true midsagittal plane x = x0, so the curve shows as it is
    ys = np.arange(curve.Y.min() - 70, curve.Y.max() + 50, h)
    zs = np.arange(curve.Z.min() - 5, curve.Z.max() + 15, h)
    Y_, Z_ = np.meshgrid(ys, zs, indexing="ij")
    pts = np.column_stack([np.full(Y_.size, curve.x0), Y_.ravel(), Z_.ravel()])
    lab_sag = label_points(pts, by_ext, cache, order).reshape(Y_.shape)
    # right: curved coronal through the centreline, spanning left-right
    ss = np.arange(curve.S[0], curve.S[-1], h)
    uu = np.arange(-45, 45, h)
    S2, U2 = np.meshgrid(ss, uu, indexing="ij")
    P = np.array([curve.point(s) for s in ss])
    pts = (P[:, None, :] + U2[..., None] * np.array([1.0, 0, 0])).reshape(-1, 3)
    lab_cor = label_points(pts, by_ext, cache, order).reshape(S2.shape)

    fig, axs = plt.subplots(1, 2, figsize=(14, 10), gridspec_kw={"width_ratios": [np.ptp(ys), 90]})
    axs[0].imshow((lab_sag + 1).T, origin="lower", cmap=cmap, vmin=0, vmax=len(order),
                  extent=(ys[0], ys[-1], zs[0], zs[-1]), aspect="equal", interpolation="nearest")
    axs[0].plot(curve.Y, curve.Z, color="w", lw=0.6, ls=":")
    cy, cz = step.centre_yz()
    for L in levels:
        c = L["vertebra_centroid_mm"]
        axs[0].text(c[1] - 30, c[2], L["vertebra_disc_naming"] + " ", color="tab:orange", ha="right",
                    va="center", fontsize=9)
    axs[0].set_xlabel("y (mm): anterior <-  -> posterior")
    axs[0].set_ylabel("z (mm)")
    axs[0].set_title("midsagittal plane x = %.1f mm (S rotates about y = %.0f, z = %.0f)" % (curve.x0, cy, cz),
                     fontsize=10)
    axs[1].imshow(lab_cor + 1, origin="lower", cmap=cmap, vmin=0, vmax=len(order),
                  extent=(uu[0], uu[-1], ss[0], ss[-1]), aspect="equal", interpolation="nearest")
    axs[1].set_xlabel("u (mm): patient right <-  -> patient left")
    axs[1].set_ylabel("arc length along the centreline (mm)")
    axs[1].set_title("curved coronal through the cord centre", fontsize=10)
    handles = [plt.Rectangle((0, 0), 1, 1, fc=cols[i + 1]) for i in range(len(order))]
    fig.legend(handles, order, loc="lower center", ncol=6, fontsize=8)
    fig.suptitle("RADO extended by %d level(s), each S = RADO's T11->T10 step (%.2f deg); orange = new "
                 "vertebrae (disc naming)" % (len(levels), step.angle_deg()), fontsize=10)
    fig.tight_layout(rect=(0, 0.05, 1, 0.97))
    fig.savefig(os.path.join(out, "extension_check.png"), dpi=110)


if __name__ == "__main__":
    main()
