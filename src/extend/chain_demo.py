"""Sympathetic chain: fused into one body per side, and a smooth-blend demo.

    src/extend/.venv/bin/python src/extend/chain_demo.py fuse
    src/extend/.venv/bin/python src/extend/chain_demo.py blend --h 0.05 --k 0.25,0.5

Reads fem/out/geometry/rado_up4 (extend_rado.py output) and writes
fem/out/geometry/rado_up4_chain_demo/.

fuse: RADO's trunk segments (neuro_sympathetic_a, straight tubes r 0.76 mm)
meet at angles and offsets, so a boolean union of them keeps a step at every
joint. The trunk is rebuilt as ONE tube along a smoothing spline through all
segment centrelines (leaving out the last 1.5 mm at each joint, where the
offsets are), then unioned with the ganglia and rami. One closed body per side.

blend: one ganglion with its two rami and the trunk through it, merged with a
smooth minimum of signed distance fields (blend radius k) and re-surfaced with
marching cubes. A boolean union keeps a sharp crease where a tube meets the
ganglion, and a crescent lip where a ramus end only partly sinks into it; the
smooth minimum fillets both. The demo box cuts the rami and trunk flat at its
faces.
"""
import argparse
import json
import os
import time

import numpy as np
from scipy.spatial import cKDTree

import meshkit as K

SRC = os.path.join(K.C.OUT, "geometry", "rado_up4")
OUT = os.path.join(K.C.OUT, "geometry", "rado_up4_chain_demo")
SPLINE_SMOOTH_MM = 0.05     # rms allowance per centreline point
JOINT_SKIP_MM = 1.5


def chain_files():
    fs = sorted(f for f in os.listdir(SRC) if f.lower().endswith(".stl")
                and K.pipe_of(os.path.splitext(f)[0]) == "sympathetic")
    xs = {f: K.read_tri(os.path.join(SRC, f)).reshape(-1, 3)[:, 0].mean() for f in fs}
    mid = 0.5 * (min(xs.values()) + max(xs.values()))
    return {"right": [f for f in fs if xs[f] < mid], "left": [f for f in fs if xs[f] >= mid]}


def trunk_path(files):
    """Smoothing spline through the trunk segments' centrelines."""
    from scipy.interpolate import splev, splprep
    trunks = sorted([f for f in files if "sympathetic_a" in f],
                    key=lambda f: K.read_tri(os.path.join(SRC, f))[:, :, 2].min())
    pts, radii = [], []
    for f in trunks:
        p, r = K.tube_centreline(K.read_tri(os.path.join(SRC, f)))
        pts.append(p)
        radii.append(r)
    keep = []
    for i, p in enumerate(pts):
        sl = np.r_[0, np.cumsum(np.linalg.norm(np.diff(p, axis=0), axis=1))]
        m = np.ones(len(p), bool)
        if i > 0:
            m &= sl > JOINT_SKIP_MM
        if i < len(pts) - 1:
            m &= sl < sl[-1] - JOINT_SKIP_MM
        keep.append(p[m])
    Q = np.vstack(keep)
    u = np.r_[0, np.cumsum(np.linalg.norm(np.diff(Q, axis=0), axis=1))]
    tck, _ = splprep(Q.T, u=u / u[-1], s=len(Q) * SPLINE_SMOOTH_MM ** 2, k=3)
    path = np.array(splev(np.linspace(0, 1, int(u[-1] / 0.25) + 1), tck)).T
    dev = cKDTree(path).query(Q)[0]
    t = np.diff(path, axis=0)
    t /= np.linalg.norm(t, axis=1)[:, None]
    bend = np.degrees(np.arccos(np.clip((t[1:] * t[:-1]).sum(1), -1, 1))) / 0.25
    return path, float(np.median(radii)), trunks, {
        "segments": len(trunks), "radius_mm": round(float(np.median(radii)), 3),
        "length_mm": round(float(np.sum(np.linalg.norm(np.diff(path, axis=0), axis=1))), 1),
        "spline_max_dev_mm": round(float(dev.max()), 3), "spline_p95_dev_mm": round(float(np.percentile(dev, 95)), 3),
        "max_bend_deg_per_mm": round(float(bend.max()), 2)}


def fuse():
    os.makedirs(OUT, exist_ok=True)
    report = {}
    for side, files in chain_files().items():
        t0 = time.time()
        path, r, trunks, info = trunk_path(files)
        np.savez(os.path.join(OUT, "trunk_%s.npz" % side), path=path, r=r)
        acc = K.manifold(*K.tube_mesh(path, r))
        others = [f for f in files if f not in trunks]
        for f in others:
            M, _ = K.load_manifold(os.path.join(SRC, f))
            acc = K.drop_slivers(acc + M)
        name = "T8-10 - neuro_sympathetic_a-chain_%s__fused.STL" % side
        assert K.pipe_of(os.path.splitext(name)[0]) == "sympathetic", name
        tri = K.manifold_tri(acc)
        K.write_stl(os.path.join(OUT, name), tri, b"fused sympathetic chain, smooth trunk")
        a = K.stlio.audit(*K.stlio.weld(K.read_tri(os.path.join(OUT, name)), tol=1e-9))
        info.update({"other_pieces": len(others), "file": name, "components": len(acc.decompose()),
                     "volume_mm3": round(acc.volume(), 1), "triangles": len(tri),
                     "open_edges": a["boundary_edges"], "nonmanifold_edges": a["nonmanifold_edges"],
                     "seconds": round(time.time() - t0, 1)})
        report[side] = info
        print(side, info, flush=True)
    json.dump(report, open(os.path.join(OUT, "fuse_report.json"), "w"), indent=1)


def surface_samples(tri, spacing):
    """Points on the surface and the triangle each lies on."""
    a, b, c = tri[:, 0], tri[:, 1], tri[:, 2]
    area = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
    cnt = np.maximum(1, np.ceil(area / spacing ** 2 * 2).astype(int))
    ti = np.repeat(np.arange(len(tri)), cnt)
    u, v = np.random.default_rng(0).random((2, len(ti)))
    flip = u + v > 1
    u[flip], v[flip] = 1 - u[flip], 1 - v[flip]
    P = np.vstack([a[ti] + u[:, None] * (b - a)[ti] + v[:, None] * (c - a)[ti], tri.reshape(-1, 3)])
    return P, np.r_[ti, np.repeat(np.arange(len(tri)), 3)]


def point_triangle_distance(p, a, b, c):
    """Exact distance from points p to triangles (a, b, c), row-wise
    (closest point by Voronoi region, Ericson, Real-Time Collision Detection 5.1.5)."""
    ab, ac, ap = b - a, c - a, p - a
    d1, d2 = (ab * ap).sum(1), (ac * ap).sum(1)
    bp = p - b
    d3, d4 = (ab * bp).sum(1), (ac * bp).sum(1)
    cp = p - c
    d5, d6 = (ab * cp).sum(1), (ac * cp).sum(1)
    va, vb, vc = d3 * d6 - d5 * d4, d5 * d2 - d1 * d6, d1 * d4 - d3 * d2
    denom = va + vb + vc
    with np.errstate(divide="ignore", invalid="ignore"):
        v = np.where(np.abs(denom) > 0, vb / denom, 0)
        w = np.where(np.abs(denom) > 0, vc / denom, 0)
        q = a + ab * v[:, None] + ac * w[:, None]                         # face interior
        t = np.clip(d1 / np.where(d1 - d3 != 0, d1 - d3, 1), 0, 1)
        e_ab = a + ab * t[:, None]
        t = np.clip(d2 / np.where(d2 - d6 != 0, d2 - d6, 1), 0, 1)
        e_ac = a + ac * t[:, None]
        t = np.clip((d4 - d3) / np.where((d4 - d3) + (d5 - d6) != 0, (d4 - d3) + (d5 - d6), 1), 0, 1)
        e_bc = b + (c - b) * t[:, None]
    q = np.where(((vc <= 0) & (d1 >= 0) & (d3 <= 0))[:, None], e_ab, q)
    q = np.where(((vb <= 0) & (d2 >= 0) & (d6 <= 0))[:, None], e_ac, q)
    q = np.where(((va <= 0) & ((d4 - d3) >= 0) & ((d5 - d6) >= 0))[:, None], e_bc, q)
    q = np.where(((d1 <= 0) & (d2 <= 0))[:, None], a, q)
    q = np.where(((d3 >= 0) & (d4 <= d3))[:, None], b, q)
    q = np.where(((d6 >= 0) & (d5 <= d6))[:, None], c, q)
    return np.linalg.norm(p - q, axis=1)


def exact_distance(pts, tri, tree, owner, k=8):
    """Exact distance to the surface: minimum over the triangles of the k
    nearest surface samples (the true nearest triangle is among them unless
    the sampling is far coarser than the triangles)."""
    _, nn = tree.query(pts, k=k, workers=-1)
    best = np.full(len(pts), np.inf)
    for j in range(k):
        t = owner[nn[:, j]]
        best = np.minimum(best, point_triangle_distance(pts, tri[t, 0], tri[t, 1], tri[t, 2]))
    return best


class Grid:
    def __init__(self, lo, hi, h):
        self.lo, self.h = lo, h
        self.shape = tuple(int(x) for x in np.ceil((hi - lo) / h).astype(int) + 1)
        self.n = int(np.prod(self.shape))

    def pts(self, idx):
        nx, ny, nz = self.shape
        i, rem = np.divmod(idx, ny * nz)
        j, k = np.divmod(rem, nz)
        return self.lo + np.column_stack([i, j, k]) * self.h


def mesh_field(tri, grid, band, chunk=2_000_000, timing=None):
    """Signed distance clamped to +-band: exact (dense samples + inside test)
    within band of the surface, sign-only beyond (from a coarse inside grid:
    a point >= band from the surface has the same side as any point within
    band - h_coarse*sqrt(3)/2 of it)."""
    t = time.time()
    P, owner = surface_samples(tri, min(0.015, grid.h / 3))
    tree = cKDTree(P)
    V, T = K.stlio.weld(tri.astype(np.float64))
    tester = K.InsideTester(V, T)
    hc = max(grid.h, band / 2)
    cg = Grid(grid.lo, grid.lo + (np.array(grid.shape) - 1) * grid.h, hc)
    inside_c = np.zeros(cg.n, bool)
    blo, bhi = P.min(0), P.max(0)
    for s in range(0, cg.n, chunk):
        idx = np.arange(s, min(s + chunk, cg.n))
        g = cg.pts(idx)
        m = ((g >= blo) & (g <= bhi)).all(1)
        if m.any():
            inside_c[idx[m]] = tester(g[m])
    t_coarse = time.time() - t
    F = np.full(grid.n, band, np.float32)
    n_band = 0
    for s in range(0, grid.n, chunk):
        idx = np.arange(s, min(s + chunk, grid.n))
        g = grid.pts(idx)
        m = ((g >= blo - band) & (g <= bhi + band)).all(1)
        if not m.any():
            continue
        gi = g[m]
        u, _ = tree.query(gi, distance_upper_bound=band, workers=-1)
        near = u < band
        n_band += int(near.sum())
        sgn = np.ones(len(gi), np.float32)
        if near.any():
            u[near] = exact_distance(gi[near], tri, tree, owner)
            sgn[np.flatnonzero(near)[tester(gi[near])]] = -1
        far = ~near
        if far.any():
            ci = np.clip(np.round((gi[far] - cg.lo) / hc).astype(int), 0, np.array(cg.shape) - 1)
            flat = (ci[:, 0] * cg.shape[1] + ci[:, 1]) * cg.shape[2] + ci[:, 2]
            sgn[far] = np.where(inside_c[flat], -1, 1)
        F[idx[m]] = sgn * np.minimum(u, band)
    if timing is not None:
        timing.append({"samples": len(P), "coarse_s": round(t_coarse, 1), "band_points": n_band,
                       "total_s": round(time.time() - t, 1)})
    return F


def tube_field(path, r, grid, band, chunk=2_000_000):
    seg = np.r_[0, np.cumsum(np.linalg.norm(np.diff(path, axis=0), axis=1))]
    ss = np.arange(0, seg[-1], min(0.01, grid.h / 5))
    dense = np.column_stack([np.interp(ss, seg, path[:, k]) for k in range(3)])
    hi = grid.lo + (np.array(grid.shape) - 1) * grid.h
    dense = dense[((dense >= grid.lo - band - r) & (dense <= hi + band + r)).all(1)]
    tree = cKDTree(dense)
    F = np.full(grid.n, band, np.float32)
    for s in range(0, grid.n, chunk):
        idx = np.arange(s, min(s + chunk, grid.n))
        u, _ = tree.query(grid.pts(idx), distance_upper_bound=band + r, workers=-1)
        F[idx] = np.clip(u - r, -band, band)
    return F


def smin(a, b, k):
    """Polynomial smooth minimum (blend width k; k = 0 is the plain union)."""
    if k <= 0:
        return np.minimum(a, b)
    hh = np.maximum(k - np.abs(a - b), 0) / k
    return np.minimum(a, b) - hh * hh * k * 0.25


def blend(h, ks, ganglion, margin):
    from skimage.measure import marching_cubes
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    side = "right" if ganglion in chain_files()["right"] else "left"
    g_tri = K.read_tri(os.path.join(SRC, ganglion))
    lo, hi = g_tri.reshape(-1, 3).min(0) - margin, g_tri.reshape(-1, 3).max(0) + margin
    pieces = {ganglion: g_tri}
    for f in chain_files()[side]:            # rami of this ganglion: non-trunk pieces touching the box
        if "sympathetic_a" in f or "ball" in f:
            continue
        tri = K.read_tri(os.path.join(SRC, f))
        P = tri.reshape(-1, 3)
        if ((P >= lo) & (P <= hi)).all(1).any():
            pieces[f] = tri
    grid = Grid(lo, hi, h)
    band = max(ks) + 2 * h
    print("box %s..%s mm, %s = %.1f M points at h = %.3f mm; band %.2f mm; pieces: %s"
          % (lo.round(1), hi.round(1), grid.shape, grid.n / 1e6, h, band, [f[8:] for f in pieces]), flush=True)
    fields, timing = [], []
    for f, tri in pieces.items():
        fields.append(mesh_field(tri, grid, band, timing=timing))
        print("  field %-42s %s" % (f[8:], timing[-1]), flush=True)
    tr = np.load(os.path.join(OUT, "trunk_%s.npz" % side))
    t = time.time()
    fields.append(tube_field(tr["path"], float(tr["r"]), grid, band))
    print("  field trunk tube %.1f s" % (time.time() - t), flush=True)
    report = {"ganglion": ganglion, "h_mm": h, "box_lo": lo.tolist(), "box_hi": hi.tolist(),
              "pieces": list(pieces), "timing": timing, "results": {}}
    for k in [0.0] + ks:
        D = fields[0].copy()
        for F in fields[1:]:
            D = smin(D, F, k)
        D = D.reshape(grid.shape)
        D[0], D[-1], D[:, 0], D[:, -1], D[:, :, 0], D[:, :, -1] = band, band, band, band, band, band
        verts, faces, _, _ = marching_cubes(D, level=0.0, spacing=(h, h, h))
        verts += lo
        a, b, c = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
        vol = np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0
        if vol < 0:
            faces, vol = faces[:, ::-1], -vol
        name = "ganglion_blend_k%.2f.stl" % k if k > 0 else "ganglion_union_k0.stl"
        tri = verts[faces]
        K.write_stl(os.path.join(OUT, name), tri, b"smooth blend demo")
        au = K.stlio.audit(*K.stlio.weld(tri, tol=1e-7))
        report["results"][name] = {"k_mm": k, "triangles": len(faces), "volume_mm3": round(float(vol), 2),
                                   "open_edges": au["boundary_edges"], "nonmanifold_edges": au["nonmanifold_edges"],
                                   "shells": au["n_shells"]}
        print("  k=%.2f: %s" % (k, report["results"][name]), flush=True)
    report["seconds"] = round(time.time() - t0, 1)
    json.dump(report, open(os.path.join(OUT, "blend_report_h%.3f.json" % h), "w"), indent=1)
    print("done in %.0f s" % (time.time() - t0))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("step", choices=["fuse", "blend"])
    ap.add_argument("--h", type=float, default=0.05, help="blend grid spacing, mm")
    ap.add_argument("--k", default="0.25,0.5", help="blend radii, mm")
    ap.add_argument("--ganglion", default="T8-10 - neuro_sympathetic_ball-1__up1.STL")
    ap.add_argument("--margin", type=float, default=3.0, help="box margin around the ganglion, mm")
    a = ap.parse_args()
    if a.step == "fuse":
        fuse()
    else:
        blend(a.h, [float(x) for x in a.k.split(",")], a.ganglion, a.margin)


if __name__ == "__main__":
    main()
