"""Extend the RADO-SCS geometry rostrally by whole vertebral levels.

    fem/.venv/bin/python src/extend/extend_rado.py                 # 4 levels, straight
    fem/.venv/bin/python src/extend/extend_rado.py --levels 3 --kyphosis-deg-per-level 2

Writes fem/out/geometry/rado_up<N>_kyph<deg>/ (git-ignored): every source STL
copied unchanged, the new STLs, manifest.json (inputs with SHA-256, every
parameter, transform and check) and extension_check.png. STL_files/ is only
read.

WHY THIS IS POSSIBLE (measured, 2026-10-09). RADO is one level unit, copied:
224 of its 232 per-level bodies are exact rigid copies of a sibling at another
level (residual ~2e-5 mm), and the 8 bodies that run the whole length (white,
grey, CSF, meninges, epidural, two blood layers, aorta) are constant
cross-sections swept along a curve. So:

  * per-level bodies: the TOP level unit (top vertebra, the disc above it, the
    roots/DRG/sleeves/vessels/ganglia at that disc) is copied N times;
  * continuous bodies: each one's perpendicular cross-section is taken just
    below its top end and swept up along an extended centreline, overlapping
    the original by --join-back mm. The original files are untouched.

Everything hangs off ONE centreline frame, fitted to the cord: copy k of the
top unit is F(p + k*D) o F(p)^-1 applied to it, where F(s) is the centreline
frame at arc length s, D the level spacing (RADO's top vertebra step, from a
fit in which the frame reproduces RADO's own vertebra-to-vertebra steps to
0.07 mm) and p the pin. The unit was shaped around RADO's curved canal, so a
copy on a straight canal matches exactly only near p; --attach roots (default)
pins it at the middle of the root complex, which keeps copied roots and
rootlets within 0.2 mm of their RADO position relative to the cord and puts
the ~1 mm residual on the vertebra, where it is electrically invisible (a gap
fills with background at the vertebra's own 0.04 S/m). Measured and recorded.

RADO's levels are rigid copies of each other but were each placed by hand
(the steps differ level to level), so the new levels are copies of RADO's TOP
level arrangement, not an average.

THE MODELLING CHOICE: --kyphosis-deg-per-level. RADO's cord tilts 17 deg from z
at the bottom and ~1 deg at the top, i.e. its kyphosis apex is at its TOP
vertebra. Repeating RADO's own top step (5.87 deg) would bend the spine
anteriorly above T10, earlier than real anatomy. The default, 0, continues
straight along the top tangent. Positive values bend in RADO's own sense.

Level names follow RADO's disc filenames (top vertebra = T10, so 4 new levels
are T9..T6), the same convention as src/freecad/README.md "Vertebral levels";
the manifest also records the paper's naming (one level higher). New discs are
named for the levels they separate (T_disk_8_9x...), other copies keep their
source name plus "__up<k>"; every new name is asserted to map to the same
tissue_map.yaml class as its source.

THE MESHER TOLERATES THE SEAMS. fem/ labels box-mesh tets by "inside any body
of this class", with one point-in-shell tester per body, so overlapping bodies
of the same class are a union and need not be conformal. What it cannot
tolerate is a gap (it becomes 0.04 S/m background); the checks look for that.
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
from scipy.optimize import minimize
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
OPEN_DIR = "open_surfaces_display_only"   # subdirectory: never scanned by config.py
RIGID_TOL_MM = 1e-3
VERT_INDEX = {**{"T%d" % i: 7 + i for i in range(1, 13)}, "L1": 20}
VERT_NAME = {v: k for k, v in VERT_INDEX.items()}


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


def rot_x_deg(R):
    """Signed rotation angle about +x of a rotation that is (nearly) pure x."""
    return float(np.degrees(np.arctan2(R[2, 1], R[1, 1])))


def apply(M, pts):
    return pts @ M[:3, :3].T + M[:3, 3]


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
    length s, (y, z) and tilt phi (tangent = (0, sin phi, cos phi)).

    Inside the cord's data range it is the cubic fit y(z); above the data it
    continues with dphi/ds = -kappa (kappa = 0: straight). RADO's own steps
    rotate +x by a positive angle going up, i.e. phi DEcreases with s, so a
    positive kappa bends in RADO's own sense."""

    def __init__(self, coef, x0, z_lo, z_hi, ext_len, kappa, ds=0.05):
        z = np.arange(z_lo, z_hi + 1e-9, ds)
        dy = np.polyval(np.polyder(coef), z)
        y = np.polyval(coef, z)
        phi = np.arctan(dy)
        seg = np.sqrt(1 + dy ** 2)
        s = np.concatenate([[0], np.cumsum(0.5 * (seg[1:] + seg[:-1]) * np.diff(z))])
        self.s_data_end = s[-1]
        se = s[-1] + np.arange(ds, ext_len + ds, ds)
        pe = phi[-1] - kappa * (se - s[-1])
        pm = 0.5 * (np.concatenate([[phi[-1]], pe[:-1]]) + pe)
        ye = y[-1] + np.cumsum(np.sin(pm) * ds)
        ze = z[-1] + np.cumsum(np.cos(pm) * ds)
        self.S = np.concatenate([s, se])
        self.Y = np.concatenate([y, ye])
        self.Z = np.concatenate([z, ze])
        self.PHI = np.concatenate([phi, pe])
        self.x0 = x0
        self._tree = cKDTree(np.column_stack([self.Y, self.Z]))

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

    def frame_map(self, sa, sb):
        """4x4 rigid map carrying the frame at sa onto the frame at sb."""
        a = self.phi(sa) - self.phi(sb)
        R = np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)]])
        M = np.eye(4)
        M[:3, :3] = R
        M[:3, 3] = self.point(sb) - R @ self.point(sa)
        return M

    def project(self, pts):
        """Arc length of the nearest centreline sample, by (y, z)."""
        pts = np.atleast_2d(pts)
        _, i = self._tree.query(pts[:, 1:3])
        return self.S[i]


# --------------------------------------------------------------------------
# Cross-section slicing, sweeping, capping
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
    """Triangulate a simple polygon (n, 2); returns index triples, or None."""
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


def cap_tris(P):
    """(triangles, method). Ear clipping; if that fails, a fan from the
    centroid, which is still correct for the ray-PARITY containment test used
    by fem/ (a fan covers each interior point an odd number of times)."""
    t = earclip(P)
    if t is not None:
        return np.array(t), "earclip"
    c = len(P)
    return np.array([(c, j, (j + 1) % c) for j in range(c)]), "fan"


def sweep_samples(curve, s0, s1, max_dphi=np.radians(0.2), max_step=10.0):
    ss = [s0]
    while ss[-1] < s1 - 1e-9:
        s = ss[-1]
        step = max_step
        while step > 0.05 and abs(curve.phi(s + step) - curve.phi(s)) > max_dphi:
            step /= 2
        ss.append(min(s + step, s1))
    return np.array(ss)


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


def clip_below(tri, p0, nrm):
    """The part of a triangle soup on the side (x - p0) . nrm < 0, with
    crossing triangles clipped exactly (RADO's swept STLs use triangles tens
    of mm long, so dropping them by centroid would leave a ragged edge)."""
    d = np.einsum("tij,j->ti", tri - p0, nrm)
    keep = tri[(d < 0).all(1)]
    out = [keep]
    for t, dd in zip(tri[~((d < 0).all(1) | (d >= 0).all(1))], d[~((d < 0).all(1) | (d >= 0).all(1))]):
        poly = []
        for i in range(3):
            a, b, da, db = t[i], t[(i + 1) % 3], dd[i], dd[(i + 1) % 3]
            if da < 0:
                poly.append(a)
            if (da < 0) != (db < 0):
                poly.append(a + da / (da - db) * (b - a))
        out.append(np.array([[poly[0], poly[j], poly[j + 1]] for j in range(1, len(poly) - 1)]).reshape(-1, 3, 3))
    return np.concatenate(out)


def cap_area_near(tri, centre, tangent, half=15.0, cos_min=0.99):
    """Area of triangles facing along the tangent within +/-half mm of
    centre: interior caps show up here, side walls do not."""
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    a = 0.5 * np.linalg.norm(n, axis=1)
    cosv = np.abs(n @ tangent) / np.maximum(2 * a, 1e-30)
    near = np.abs((tri.mean(1) - centre) @ tangent) < half
    return float(a[(cosv > cos_min) & near].sum())


def curve_frames(curve, s_samples):
    """(centre, e_u, e_v) along the centreline: e_u = +x, e_v = in-plane normal."""
    return [(curve.point(s), np.array([1.0, 0, 0]), curve.normal(s)) for s in s_samples]


def line_frames(p0, d, eu, ev, lengths):
    return [(p0 + t * d, eu, ev) for t in lengths]


def sweep_body(frames, loops2d, caps=("bottom", "top")):
    """Closed tube per loop (side walls + two caps), all loops in one array.
    Nested loops (holes) are handled by the parity containment test. Pass
    caps=("top",) for the open display/size-field surface."""
    tris, methods = [], []
    for P in loops2d:
        n = len(P)
        ring = np.array([[c + P[j, 0] * eu + P[j, 1] * ev
                          for j in range(n)] for c, eu, ev in frames])   # (ns, n, 3)
        for i in range(len(frames) - 1):
            for j in range(n):
                k = (j + 1) % n
                tris.append([ring[i, j], ring[i, k], ring[i + 1, k]])
                tris.append([ring[i, j], ring[i + 1, k], ring[i + 1, j]])
        ct, how = cap_tris(P)
        methods.append(how)
        ends = [e for e, name in ((0, "bottom"), (len(frames) - 1, "top")) if name in caps]
        for end in ends:
            pts = ring[end]
            if how == "fan":
                pts = np.vstack([pts, pts.mean(0)])
            for a, b, c in ct:
                tris.append([pts[a], pts[b], pts[c]] if end else [pts[a], pts[c], pts[b]])
    return np.array(tris), methods


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


def section_grid(curve, s, half_u, v_lo, v_hi, h):
    u = np.arange(-half_u, half_u + h / 2, h)
    v = np.arange(v_lo, v_hi + h / 2, h)
    U, W = np.meshgrid(u, v, indexing="ij")
    pts = curve.point(s) + U.ravel()[:, None] * np.array([1.0, 0, 0]) + W.ravel()[:, None] * curve.normal(s)
    return pts, U.shape


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


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--levels", type=int, default=4, help="vertebral levels to add above the top")
    ap.add_argument("--kyphosis-deg-per-level", type=float, default=0.0,
                    help="bending per added level, in RADO's own sense (0 = straight; "
                         "RADO's own top step is ~5.87)")
    ap.add_argument("--spacing", default="top", help="level spacing: 'top' (RADO's top step), "
                    "'mean' (mean of RADO's steps) or a number in mm of arc length")
    ap.add_argument("--join-back", type=float, default=10.0,
                    help="mm below each continuous body's top end where its section is "
                         "taken and its extension starts (= overlap with the original). "
                         "RADO's continuous bodies are not the constant section within "
                         "~3-8 mm of their ends (measured), so keep this >= 10")
    ap.add_argument("--attach", choices=["roots", "vertebra", "fit"], default="roots",
                    help="where the copied top unit is pinned to the centreline (see docstring)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--no-check", action="store_true", help="skip the section checks and figure")
    args = ap.parse_args()
    t0 = time.time()
    K = args.levels
    warnings = []
    out = args.out or os.path.join(C.OUT, "geometry", "rado_up%d_kyph%g" % (K, args.kyphosis_deg_per_level))
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

    # ---- centreline ----------------------------------------------------------
    cpaths = C.TISSUE_BODIES["white"] + C.TISSUE_BODIES["grey"]
    cen = cord_centroids(cpaths)
    coef = np.polyfit(cen[:, 2], cen[:, 1], 3)
    fit_res = cen[:, 1] - np.polyval(coef, cen[:, 2])
    x0 = float(cen[:, 0].mean())
    ext_len = K * 40.0 + 120.0
    curve = Curve(coef, x0, cen[0, 2] - 12.0, cen[-1, 2], ext_len, 0.0)
    print("centreline: cubic y(z) residual rms %.4f mm max %.4f mm; x drift sd %.4f mm"
          % (np.sqrt(np.mean(fit_res ** 2)), np.abs(fit_res).max(), cen[:, 0].std()))

    # ---- vertebra attachments: make the frame reproduce RADO's steps ------
    verts = sorted([b for b in per_level if b["pipe"] == "vertebra"], key=lambda b: b["centroid"][2])
    P = [b["tri"].reshape(-1, 3) for b in verts]
    steps_measured = []
    for a, b in zip(P[:-1], P[1:]):
        R, t, res = rigid_fit(a, b)
        assert res < RIGID_TOL_MM, "vertebrae are no longer exact rigid copies"
        steps_measured.append({"rot_x_deg": rot_x_deg(R), "R": R, "t": t})
    sub = np.arange(0, len(P[0]), max(1, len(P[0]) // 3000))
    Ps = [p[sub] for p in P]

    def objective(sv):
        e = 0.0
        for i in range(len(Ps) - 1):
            e += np.mean(np.sum((apply(curve.frame_map(sv[i], sv[i + 1]), Ps[i]) - Ps[i + 1]) ** 2, 1))
        return e
    s_init = np.array([curve.project(p.mean(0))[0] for p in P])
    fit = minimize(objective, s_init, method="Nelder-Mead",
                   options={"xatol": 1e-5, "fatol": 1e-12, "maxiter": 20000})
    s_att = fit.x
    fit_steps = []
    for i in range(len(P) - 1):
        M = curve.frame_map(s_att[i], s_att[i + 1])
        r = np.linalg.norm(apply(M, P[i]) - P[i + 1], axis=1)
        fit_steps.append({"from": verts[i]["stem"], "to": verts[i + 1]["stem"],
                          "rado_rot_x_deg": steps_measured[i]["rot_x_deg"],
                          "frame_rot_x_deg": float(np.degrees(curve.phi(s_att[i]) - curve.phi(s_att[i + 1]))),
                          "arc_length_mm": float(s_att[i + 1] - s_att[i]),
                          "residual_rms_mm": float(np.sqrt(np.mean(r ** 2))),
                          "residual_max_mm": float(r.max())})
        print("  frame vs RADO step %s->%s: RADO %.3f deg, frame %.3f deg, arc %.2f mm, residual rms %.3f max %.3f mm"
              % (verts[i]["stem"][8:], verts[i + 1]["stem"][8:], fit_steps[-1]["rado_rot_x_deg"],
                 fit_steps[-1]["frame_rot_x_deg"], fit_steps[-1]["arc_length_mm"],
                 fit_steps[-1]["residual_rms_mm"], fit_steps[-1]["residual_max_mm"]))
        if fit_steps[-1]["residual_max_mm"] > 0.5:
            warnings.append("frame reproduces RADO step %d only to %.2f mm max" % (i, fit_steps[-1]["residual_max_mm"]))
    s3 = float(s_att[-1])
    if args.spacing == "top":
        D = float(s_att[-1] - s_att[-2])
    elif args.spacing == "mean":
        D = float(np.mean(np.diff(s_att)))
    else:
        D = float(args.spacing)

    # rebuild the curve with the requested extension curvature
    kappa = np.radians(args.kyphosis_deg_per_level) / D
    curve = Curve(coef, x0, cen[0, 2] - 12.0, cen[-1, 2], ext_len, kappa)

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
    # CENTROID (the fitted attachment s3 is a frame position, not a centre).
    # RADO's levels are rigid copies but were NOT placed by one repeated
    # transform (e.g. the left DRG steps 31.3 / 27.6 / 25.8 mm with sideways
    # shifts while the vertebrae step 31.9 / 30.8 mm), so the TOP level is
    # copied as one rigid snapshot of RADO's own top-level arrangement. Each
    # group's own top step is recorded for reference, not as a pass/fail.
    s_vref = float(curve.project(verts[-1]["centroid"])[0])
    copy_set, group_report = [], []
    for key, mem in sorted(groups.items()):
        sc = np.array([curve.project(b["centroid"])[0] for b in mem])
        top = sc.max()
        offset = (top - s_vref) / D
        tops = [b for b, s in zip(mem, sc) if s > top - 0.5 * D]
        include = offset > -0.5
        if -0.65 < offset < -0.35:
            warnings.append("group %s top member is ambiguously placed (%.2f levels from the top vertebra)" % (key, offset))
        own_step = None
        lower = [b for b, s in zip(mem, sc) if top - 1.5 * D <= s < top - 0.5 * D]
        if lower:
            tb = tops[0]
            cand = [lb for lb in lower if lb["tri"].shape == tb["tri"].shape]
            if cand:
                lb = min(cand, key=lambda b: abs(b["centroid"][0] - tb["centroid"][0]))   # same side
                R, t, r = rigid_fit(lb["tri"].reshape(-1, 3), tb["tri"].reshape(-1, 3))
                if r < RIGID_TOL_MM:
                    own_step = {"rot_x_deg": round(rot_x_deg(R), 3), "translation_mm": t.round(3).tolist()}
        group_report.append({"group": key, "class": mem[0]["pipe"], "members": len(mem),
                             "top_members": [b["stem"] for b in tops],
                             "top_offset_levels_from_top_vertebra_centroid": round(float(offset), 3),
                             "copied": bool(include), "rado_own_top_step": own_step})
        if include:
            copy_set.extend(tops)
    print("per-level groups: %d (%d copied as the top-level unit, %d bodies per new level)"
          % (len(groups), sum(g["copied"] for g in group_report), len(copy_set)))

    # ---- where the copied unit is pinned ------------------------------------
    # RADO's top unit was shaped around a CURVED canal; above RADO the canal is
    # straight (or differently curved), so the copy cannot match everywhere.
    # Wherever it is pinned, it matches exactly there and drifts by roughly
    # curvature x distance^2 / 2 away from it. Pinned at the fitted vertebra
    # attachment, rootlets ended up to 3.5 mm off their RADO position relative
    # to the cord; pinned at the middle of the root complex, 0.2 mm, with the
    # vertebra ~1 mm off instead. Bone is the right place to put the error:
    # any gap it opens is filled with background at 0.04 S/m, the vertebra's
    # own conductivity, while roots, rootlets and their CSF sleeves are
    # stimulation targets. Both offsets are measured and recorded below.
    root_cx = [b for b in copy_set if b["pipe"] in ("root", "drg", "csf", "dura")]
    rc_pts = np.vstack([b["tri"].reshape(-1, 3)[::7] for b in root_cx])
    s_roots = float(np.median(curve.project(rc_pts)))
    s_attach = {"roots": s_roots, "vertebra": s_vref, "fit": s3}[args.attach]

    def radial_min(P):
        sp = curve.project(P)
        c = np.array([curve.point(x) for x in sp])
        n = np.array([curve.normal(x) for x in sp])
        r = P - c
        return np.hypot(r[:, 0], np.einsum("ij,ij->i", r, n))
    attach_check = {}
    M2 = curve.frame_map(s_attach, s_attach + 2 * D)
    for label, cls in (("roots_and_rootlets", ("root",)), ("drg", ("drg",)), ("vertebra", ("vertebra",))):
        bs = [b for b in copy_set if b["pipe"] in cls]
        d = [float(radial_min(apply(M2, b["tri"].reshape(-1, 3))).min() - radial_min(b["tri"].reshape(-1, 3)).min())
             for b in bs]
        attach_check[label] = {"bodies": len(bs), "innermost_radius_change_mean_mm": round(float(np.mean(d)), 3),
                               "innermost_radius_change_max_abs_mm": round(float(np.max(np.abs(d))), 3)}
    print("pinned at '%s' (s = %.2f; roots %.2f, vertebra centroid %.2f, fit %.2f). Copy vs RADO, distance "
          "to the cord axis: roots/rootlets max %.2f mm, DRG max %.2f mm, vertebra max %.2f mm"
          % (args.attach, s_attach, s_roots, s_vref, s3,
             attach_check["roots_and_rootlets"]["innermost_radius_change_max_abs_mm"],
             attach_check["drg"]["innermost_radius_change_max_abs_mm"],
             attach_check["vertebra"]["innermost_radius_change_max_abs_mm"]))
    if attach_check["roots_and_rootlets"]["innermost_radius_change_max_abs_mm"] > 0.5:
        warnings.append("copied roots/rootlets sit up to %.2f mm off their RADO position relative to the cord"
                        % attach_check["roots_and_rootlets"]["innermost_radius_change_max_abs_mm"])

    # ---- level names -------------------------------------------------------
    top_disc = max(groups["discs"], key=lambda b: b["centroid"][2])
    m = re.search(r"T_disk_(\d+)_(\d+|L)", top_disc["stem"])
    top_vert_disc_naming = int(m.group(2))          # the vertebra below the top disc
    levels = []
    for k in range(1, K + 1):
        v = top_vert_disc_naming - k
        levels.append({"k": k, "vertebra_disc_naming": "T%d" % v, "vertebra_paper_naming": "T%d" % (v - 1),
                       "disc_above_disc_naming": "T%d/T%d" % (v - 1, v),
                       "pin_arc_length_mm": s_attach + k * D,
                       "vertebra_centroid_arc_length_mm": float(curve.project(apply(
                           curve.frame_map(s_attach, s_attach + k * D), verts[-1]["centroid"][None])[0])[0]),
                       "transform_from_top_level": curve.frame_map(s_attach, s_attach + k * D).round(9).tolist(),
                       "rot_x_deg_from_top_level": float(np.degrees(curve.phi(s_attach) - curve.phi(s_attach + k * D)))})

    # ---- copy everything, write the new per-level bodies -------------------
    outputs = []
    for b in bodies:
        shutil.copy2(b["path"], os.path.join(out, b["file"]))
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
            tri = apply(M, b["tri"].reshape(-1, 3)).reshape(-1, 3, 3)
            path = os.path.join(out, stem + ".STL")
            write_stl(path, tri)
            outputs.append({"file": stem + ".STL", "kind": "copy", "level_k": k, "source": b["file"],
                            "class": b["pipe"], "triangles": len(tri)})
    print("wrote %d level copies" % len(outputs))

    # ---- sweep the continuous bodies --------------------------------------
    # Canal contents are constant sections in the centreline frame, so they are
    # swept along the extended centreline. A body that is NOT constant in that
    # frame (the aorta runs on its own path) is extruded straight along its
    # own top axis instead. The choice is made by measurement, per body.
    sweep_report, open_outputs = [], []
    dz_new = float(curve.point(curve.s_data_end + K * D)[2] - curve.point(curve.s_data_end)[2])
    for b in continuous:
        V, T = stlio.weld(b["tri"])
        s_top = float(curve.project(V).max())
        s_join = s_top - args.join_back
        s_end = s_top + K * D
        ex = np.array([1.0, 0, 0])
        loops_hi, dropped = section_loops(V, T, curve.point(s_join), ex, curve.normal(s_join), curve.tangent(s_join))
        loops_lo, _ = section_loops(V, T, curve.point(s_join - 20), ex, curve.normal(s_join - 20),
                                    curve.tangent(s_join - 20))
        constancy = loops_dice(loops_hi, loops_lo)
        drift = float(np.linalg.norm(loops_centroid(loops_hi) - loops_centroid(loops_lo)))
        if drift < 1.0:
            method = "centreline"
            loops2 = loops_hi
            samples = sweep_samples(curve, s_join, s_end)
            frames = curve_frames(curve, samples)
            sm = s_join + 0.5 * args.join_back
            mid = (curve.point(sm), ex, curve.normal(sm), curve.tangent(sm))
        else:
            method = "own_axis"
            z_top = float(V[:, 2].max())
            p_axis, d = own_axis(V, T, z_top)
            p0 = p_axis + d * ((z_top - args.join_back - p_axis[2]) / d[2])
            eu = ex - (ex @ d) * d
            eu /= np.linalg.norm(eu)
            ev = np.cross(d, eu)
            loops2, dropped = section_loops(V, T, p0, eu, ev, d)
            length = (z_top + dz_new - p0[2]) / d[2]
            samples = np.array([0.0, length])
            frames = line_frames(p0, d, eu, ev, samples)
            mid = (p0 + 0.5 * args.join_back / d[2] * d, eu, ev, d)
        tri, methods = sweep_body(frames, loops2)
        stem = "%s__sweep_up%d" % (b["stem"], K)
        assert class_of(stem, entries) == b["cls"], stem
        write_stl(os.path.join(out, stem + ".STL"), tri)
        # Open DISPLAY / SIZE-FIELD surface: the original cut at the join
        # plane + the sweep's walls and top end. No interior caps, no doubled
        # walls in the overlap. Written to a SUBDIRECTORY on purpose:
        # config.py classifies every .stl at the top level of the geometry
        # directory, and an open surface there would become a tissue body and
        # break the parity containment test. Never classify these.
        c_j, n_j = frames[0][0], np.cross(frames[0][1], frames[0][2])
        n_j = n_j / np.linalg.norm(n_j)
        if n_j @ (frames[-1][0] - frames[0][0]) < 0:
            n_j = -n_j
        walls_top, _ = sweep_body(frames, loops2, caps=("top",))
        open_tri = np.concatenate([clip_below(b["tri"], c_j, n_j), walls_top])
        os.makedirs(os.path.join(out, OPEN_DIR), exist_ok=True)
        open_name = "%s__open_surface_up%d.STL" % (b["stem"], K)
        write_stl(os.path.join(out, OPEN_DIR, open_name), open_tri)
        cap_closed = (cap_area_near(b["tri"], c_j, n_j) + cap_area_near(tri, c_j, n_j))
        cap_open = cap_area_near(open_tri, c_j, n_j)
        sweep_report.append({"source": b["file"], "class": b["pipe"], "file": stem + ".STL",
                             "method": method,
                             "constancy_dice_in_centreline_frame_20mm": round(constancy, 4),
                             "section_centroid_drift_in_centreline_frame_20mm": round(drift, 4),
                             "loops": len(loops2), "loop_points": [len(q) for q in loops2],
                             "dropped_tiny_loops": dropped, "caps": methods,
                             "s_join_mm": s_join, "s_end_mm": s_end, "overlap_mm": args.join_back,
                             "path_samples": len(samples), "triangles": len(tri),
                             "mid_overlap_frame": {"centre": mid[0].tolist(), "e_u": mid[1].tolist(),
                                                   "e_v": mid[2].tolist(), "normal": mid[3].tolist()},
                             "open_surface": OPEN_DIR + "/" + open_name,
                             "cap_area_near_join_mm2": {"closed_pair": round(cap_closed, 3),
                                                        "open_surface": round(cap_open, 3)},
                             "loops_uv": [q.round(5).tolist() for q in loops2]
                             if b["pipe"] == "epidural" else None})
        outputs.append({"file": stem + ".STL", "kind": "sweep", "source": b["file"], "class": b["pipe"],
                        "triangles": len(tri)})
        open_outputs.append({"file": OPEN_DIR + "/" + open_name, "kind": "open_surface_display_size_only",
                             "source": b["file"], "class": b["pipe"], "triangles": len(open_tri)})
        if cap_open > 0.01 * max(cap_closed, 1e-9) + 0.05:
            warnings.append("%s: open surface still has %.2f mm2 of cap-like area near the join" % (open_name, cap_open))
        if "fan" in methods:
            warnings.append("%s: a cap used the fan fallback (parity-correct, geometrically folded)" % stem)
        print("  %-10s %-28s section over 20 mm: Dice %.3f, centroid drift %.3f mm; %d loop(s), caps %s"
              % (method, b["stem"][8:], constancy, drift, len(loops2), "/".join(sorted(set(methods)))))

    for o in outputs + open_outputs:
        o["sha256"] = sha256(os.path.join(out, o["file"]))
    print("open display/size-field surfaces in %s/: interior cap area near the joins %.1f mm2 "
          "(closed pairs) -> %.2f mm2 (open)" % (OPEN_DIR,
          sum(r["cap_area_near_join_mm2"]["closed_pair"] for r in sweep_report),
          sum(r["cap_area_near_join_mm2"]["open_surface"] for r in sweep_report)))

    # ---- checks ------------------------------------------------------------
    checks = {}
    if not args.no_check:
        checks = run_checks(curve, out, bodies, outputs, levels, s_vref, D, sweep_report, warnings)

    manifest = {
        "tool": "src/extend/extend_rado.py",
        "created": time.strftime("%Y-%m-%d %H:%M:%S"),
        "parameters": {"levels": K, "kyphosis_deg_per_level": args.kyphosis_deg_per_level,
                       "spacing": args.spacing, "join_back_mm": args.join_back},
        "inputs": {b["file"]: sha256(b["path"]) for b in bodies},
        "tissue_map_sha256": sha256(C.TISSUE_MAP_YAML),
        "centreline": {"model": "x = x0; y = cubic(z) inside the cord, then dphi/ds = -kappa",
                       "x0_mm": x0, "cubic_coef_high_first": coef.tolist(),
                       "fit_residual_rms_mm": float(np.sqrt(np.mean(fit_res ** 2))),
                       "fit_residual_max_mm": float(np.abs(fit_res).max()),
                       "x_drift_sd_mm": float(cen[:, 0].std()),
                       "data_z_range_mm": [float(cen[0, 2]), float(cen[-1, 2])],
                       "tilt_deg_at_data_top": float(np.degrees(curve.phi(curve.s_data_end))),
                       "kappa_rad_per_mm": float(kappa)},
        "vertebra_attachment_fit": {"arc_lengths_mm": s_att.tolist(), "steps": fit_steps,
                                    "objective": float(fit.fun)},
        "level_spacing_arc_mm": D,
        "pin": {"mode": args.attach, "arc_length_mm": s_attach, "roots_mm": s_roots,
                "vertebra_centroid_mm": s_vref, "vertebra_fit_mm": s3, "check_level_k2": attach_check},
        "new_levels": levels,
        "level_naming_note": "disc naming (default, used in new disc filenames) makes RADO's "
                             "vertebrae T10-T12; the paper's naming is one level higher. "
                             "See src/freecad/README.md 'Vertebral levels'.",
        "per_level_groups": group_report,
        "sweeps": sweep_report,
        "outputs": outputs,
        "open_surfaces_display_and_size_field_only": {
            "dir": OPEN_DIR,
            "note": "NOT closed, NEVER classify. One surface per continuous tissue: the "
                    "original cut at the join plane + the sweep's walls and top end. For "
                    "visualiser surfaces and mesh size-field clouds; the closed bodies at "
                    "the top level remain the classification set.",
            "files": open_outputs},
        "checks": checks,
        "warnings": warnings,
        "elapsed_s": round(time.time() - t0, 1),
    }
    with open(os.path.join(out, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1)
    print("\n%d warning(s)" % len(warnings))
    for w in warnings:
        print("  - " + w)
    print("wrote %s (%d new STLs) in %.0f s" % (out, len(outputs), time.time() - t0))


def run_checks(curve, out, bodies, outputs, levels, s_vref, D, sweep_report, warnings):
    """Three checks, each aimed at one way the extension could be wrong.

    1. seams: the swept extension against the original body, mid-overlap.
    2. copies: RADO's top-level unit + its continuous bodies against copy k +
       the sweeps, on the same frame-relative sections. Only the copied bodies
       take part: RADO's levels were each placed individually, so comparing
       whole levels compares different arrangements, not copy fidelity.
    3. gaps: in the new region, every point inside the canal's outer outline
       must belong to some tissue (a gap would become 0.04 S/m background).
    """
    entries = parse_tissue_map(C.TISSUE_MAP_YAML)
    order = [p for p, _ in C.ANATOMY_CLASS]
    cache = TesterCache()

    def index(files):
        by = {}
        for fn in files:
            path = os.path.join(out, fn)
            pipe = PIPE_OF_CLASS.get(class_of(os.path.splitext(fn)[0], entries))
            if not pipe:
                continue
            v = stlio.read_stl(path)[0]
            by.setdefault(pipe, []).append((path, (v.min(0), v.max(0))))
        return by
    res = {}

    # 1. seams ---------------------------------------------------------------
    seams = []
    for sw in sweep_report:
        f = sw["mid_overlap_frame"]
        c, eu, ev, nrm = (np.array(f[k]) for k in ("centre", "e_u", "e_v", "normal"))
        src, ext = os.path.join(out, sw["source"]), os.path.join(out, sw["file"])
        V, T = stlio.read_stl(src)
        loops, _ = section_loops(V, T, c, eu, ev, nrm)
        allp = np.vstack(loops)
        h = 0.05
        u = np.arange(allp[:, 0].min() - .5, allp[:, 0].max() + .5, h)
        w = np.arange(allp[:, 1].min() - .5, allp[:, 1].max() + .5, h)
        U, W = np.meshgrid(u, w, indexing="ij")
        pts = c + U.ravel()[:, None] * eu + W.ravel()[:, None] * ev
        a, b = cache.get(src)(pts), cache.get(ext)(pts)
        seams.append({"body": sw["source"], "method": sw["method"], "dice": round(dice(a, b), 4),
                      "area_original_mm2": round(float(a.sum() * h * h), 3),
                      "area_extension_mm2": round(float(b.sum() * h * h), 3),
                      "xor_area_mm2": round(float((a ^ b).sum() * h * h), 4)})
    res["seams"] = seams
    print("seams (Dice, mismatched area): " + ", ".join(
        "%s %.3f/%.3fmm2" % (s["body"][8:].split("-")[0][:13], s["dice"], s["xor_area_mm2"]) for s in seams))

    # 2. copies --------------------------------------------------------------
    # Label a section through RADO's top unit (+ the original continuous
    # bodies), then label the SAME points carried by M_k in copy k + sweeps.
    # Copies are exact, so disagreement can only come from the sweeps (canal
    # straight above RADO, curved inside it) or from a transform error.
    cont = [sw["source"] for sw in sweep_report]
    sweeps = [sw["file"] for sw in sweep_report]
    unit_src = sorted({o["source"] for o in outputs if o["kind"] == "copy"})
    by_ref = index(unit_src + cont)
    s_join_min = min(sw["s_join_mm"] for sw in sweep_report)
    reps = []
    for L in levels:
        M = np.array(L["transform_from_top_level"])
        by_k = index([o["file"] for o in outputs if o["kind"] == "copy" and o["level_k"] == L["k"]] + sweeps)
        for f in (-0.25, 0.0, 0.25, 0.5, 0.75):
            p0, _ = section_grid(curve, s_vref + f * D, 45, -45, 40, 0.2)
            l0 = label_points(p0, by_ref, cache, order)
            l1 = label_points(apply(M, p0), by_k, cache, order)
            occ = (l0 >= 0) | (l1 >= 0)
            per_level_cls = np.isin(l0, [order.index(c) for c in
                                         ("drg", "root", "sympathetic", "disc", "vertebra")])
            reps.append({"level_k": L["k"], "offset_levels": f,
                         "agreement_where_occupied": round(float((l0 == l1)[occ].mean()), 5),
                         "agreement_on_per_level_tissue": round(float((l0 == l1)[per_level_cls].mean()), 5)
                         if per_level_cls.any() else None,
                         "points_occupied": int(occ.sum())})
        # Gate on per-level tissue only: continuous bodies follow the LOCAL
        # frame, so carrying them with the top unit's M_k shifts them by
        # (distance from centreline) x (RADO's curvature between the section
        # and s3); the aorta, 40 mm anterior, moves most. Sections beyond the
        # sweep join compare against RADO's own non-constant end zone.
        mine = [r for r in reps if r["level_k"] == L["k"]]
        for r in mine:
            r["reference_in_rado_end_zone"] = bool(s_vref + r["offset_levels"] * D > s_join_min)
        worst_pl = min(r["agreement_on_per_level_tissue"] for r in mine if r["agreement_on_per_level_tissue"] is not None)
        worst = min(r["agreement_where_occupied"] for r in mine if not r["reference_in_rado_end_zone"])
        print("  copy check level +%d (%s): per-level tissue agreement >= %.4f; all tissue >= %.4f "
              "(outside RADO's end zone)" % (L["k"], L["vertebra_disc_naming"], worst_pl, worst))
        if worst_pl < 0.98:
            warnings.append("copy check level +%d: per-level tissue agreement %.4f" % (L["k"], worst_pl))
    res["copies"] = reps

    # 3. gaps ----------------------------------------------------------------
    epi = next(sw for sw in sweep_report if sw["class"] == "epidural")
    loops = [np.array(q) for q in epi["loops_uv"]]
    outer = max(loops, key=lambda q: abs(area2(q)))
    by_ext = index([b["file"] for b in bodies] + [o["file"] for o in outputs])
    gaps = []
    s_hi = s_vref + (len(levels) + 0.75) * D
    s_join_min = min(sw["s_join_mm"] for sw in sweep_report)
    base = list(np.arange(s_vref - 2.0 * D, s_join_min - 2.0, 5.0))   # RADO's own canal
    for s in base + list(np.arange(s_join_min - 2.0, s_hi, 5.0)):
        lo, hi = outer.min(0), outer.max(0)
        U, W = np.meshgrid(np.arange(lo[0], hi[0], 0.1), np.arange(lo[1], hi[1], 0.1), indexing="ij")
        inside = loops_mask([outer], U, W).ravel()
        pts = curve.point(s) + U.ravel()[inside, None] * np.array([1.0, 0, 0]) + W.ravel()[inside, None] * curve.normal(s)
        lab = label_points(pts, by_ext, cache, order)
        gaps.append({"s_mm": round(float(s), 2), "z_mm": round(float(curve.point(s)[2]), 2),
                     "region": "rado" if s < s_join_min - 2.0 else "extension",
                     "canal_points": int(inside.sum()), "background_points": int((lab == -1).sum())})
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
        plot(curve, by_ext, cache, order, entries, levels, s_vref, D, out)
        res["figure"] = "extension_check.png"
    except Exception as e:      # the figure is a convenience, never a gate
        warnings.append("figure failed: %r" % (e,))
    return res


def plot(curve, by_ext, cache, order, entries, levels, s3, D, out):
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
    h = 0.3
    s_lo, s_hi = curve.S[0], s3 + (len(levels) + 0.9) * D
    ss = np.arange(s_lo, s_hi, h)
    # midsagittal plane x = x0, in (s, v) coordinates
    vv = np.arange(-50, 35, h)
    S_, V_ = np.meshgrid(ss, vv, indexing="ij")
    P = np.array([curve.point(s) for s in ss])
    N = np.array([curve.normal(s) for s in ss])
    pts = (P[:, None, :] + V_[..., None] * N[:, None, :]).reshape(-1, 3)
    lab_sag = label_points(pts, by_ext, cache, order).reshape(S_.shape)
    # curved coronal: through the centreline, spanning left-right
    uu = np.arange(-45, 45, h)
    S2, U2 = np.meshgrid(ss, uu, indexing="ij")
    pts = (P[:, None, :] + U2[..., None] * np.array([1.0, 0, 0])).reshape(-1, 3)
    lab_cor = label_points(pts, by_ext, cache, order).reshape(S2.shape)

    fig, axs = plt.subplots(1, 2, figsize=(13, 9), gridspec_kw={"width_ratios": [85, 90]})
    for ax, lab, x_ext, xl in ((axs[0], lab_sag, (vv[0], vv[-1]), "v (mm): anterior <-  -> posterior"),
                               (axs[1], lab_cor, (uu[0], uu[-1]), "u (mm): patient right <-  -> patient left")):
        ax.imshow(lab + 1, origin="lower", cmap=cmap, vmin=0, vmax=len(order),
                  extent=(x_ext[0], x_ext[1], ss[0], ss[-1]), aspect="equal", interpolation="nearest")
        ax.axhline(curve.s_data_end, color="w", ls=":", lw=0.8)
        for L in levels:
            ax.text(x_ext[1], L["vertebra_centroid_arc_length_mm"], " " + L["vertebra_disc_naming"], color="tab:orange",
                    va="center", fontsize=8)
        ax.text(x_ext[1], s3, " T%s (RADO top)" % (int(levels[0]["vertebra_disc_naming"][1:]) + 1),
                color="w", va="center", fontsize=8)
        ax.set_xlabel(xl)
        ax.set_ylabel("arc length along centreline (mm)")
    axs[0].set_title("midsagittal (straightened along the centreline)", fontsize=10)
    axs[1].set_title("curved coronal through the cord centre", fontsize=10)
    handles = [plt.Rectangle((0, 0), 1, 1, fc=cols[i + 1]) for i in range(len(order))]
    fig.legend(handles, order, loc="lower center", ncol=6, fontsize=8)
    fig.suptitle("RADO extended by %d level(s); dotted line = top of RADO's cord data; "
                 "orange = new vertebral levels (disc naming)" % len(levels), fontsize=10)
    fig.tight_layout(rect=(0, 0.05, 1, 0.97))
    fig.savefig(os.path.join(out, "extension_check.png"), dpi=110)


if __name__ == "__main__":
    main()
