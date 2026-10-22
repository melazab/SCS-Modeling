"""Continuous 3D tract bodies that tile the extended model's cord.

    python src/atlas/build_tract_bodies.py [--atlas DIR] [--geometry DIR] [--out DIR]

Needs numpy, scipy, scikit-image, manifold3d; reads fem/scripts/inside.py and stlio.py
(imported, not modified). Runs on Pioneer in ~/scs/atlas-venv (~/scs/tract_bodies.sbatch).

register_tracts.py gives one 2D tract map per vertebra (T6-T12), all on RADO's
constant cord section in the centreline frame (u = patient left, v =
posterior). This builds one closed surface per tract and side:

  1. Grid: (u, v) every GRID_MM over the cord's section, s every DS_MM along
     FIBRE_CENTRELINE. The cord is the extended model's own white + grey
     matter STLs (inside test at every grid point), so the bodies start and
     end where the anatomy's cord does, tilted end planes included. (The first
     build, 2026-10-21, filled the whole centreline and ran ~10 mm past the
     cord's caudal end; Mohamed saw it.)
  2. Tract probabilities are interpolated linearly in arc length between the
     two nearest vertebrae (constant beyond T6 and T12); each cord voxel takes
     the most probable tract. Cord voxels outside RADO's sampled section (a
     voxel's worth at the edge) take the nearest labelled voxel's tract.
  3. Islands: a 6-connected piece of a tract smaller than ISLAND_FRACTION of
     that tract's largest piece is relabelled to its nearest other tract.
     Corner-only contacts are what made non-manifold edges in the first build.
  4. Shared surfaces: every label (tracts + outside) gets a smoothed indicator
     s_k (Gaussian, SMOOTH_VOX); tract i's surface is the zero set of
     s_i - max_{j != i} s_j. Two neighbours compute the same function with
     opposite sign, so marching cubes puts their shared boundary in the same
     place: no gaps and no overlaps between bodies (the first build smoothed
     each tract alone and left them 5 % short of the cord's area).
  5. Vertices go to 3D as c(s) + u U(s) + v V(s): the frame the maps were
     sampled in, so the bodies sit in the extended model's cord.
  6. Ends: the cord's ends are planes oblique to its axis (caudal z = 59.43,
     20 deg; cranial 18.6 deg), and a 0.5 mm grid turned them into terraces
     (Mohamed, 2026-10-22: "why is the other end blocky looking?"). So the
     grid runs END_EXT_MM past both ends with the full cord section carried
     out, and each body is trimmed by the two end planes fitted to the cord
     STLs' own end caps (manifold3d trim_by_plane): flat ends, flush with the
     cord's.

They are visual and fiber-placement aids only; the FEM never reads them.
Writes <out>/<id>_<side>_<short>.stl (binary, mm) and <out>/bodies.json.
"""
import argparse
import glob
import json
import os
import struct
import sys
import time

import manifold3d as m3d
import numpy as np
from scipy.ndimage import distance_transform_edt, gaussian_filter, label as cc_label
from skimage.measure import marching_cubes

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "fem", "scripts"))
sys.path.insert(0, os.path.join(ROOT, "src", "freecad"))
from inside import InsideTester  # noqa: E402
from stlio import read_stl  # noqa: E402
import scs_levels  # noqa: E402

LEVELS = ("T12", "T11", "T10", "T9", "T8", "T7", "T6")
# Each level's map sits at its vertebra's volume centroid projected on the
# centreline (src/freecad/scs_levels.py, from the STLs). The first builds used
# the vertebra labels' z, up to 8.5 mm off (labels mark transverse-process tips).
DS_MM = 0.5                    # along the cord (PAM50's slice spacing)
GRID_MM = 0.1                  # in-plane; the maps are 0.05 mm, a multiple of this
U_HALF, V_HALF = 4.8, 3.2      # mm; RADO's cord is 8.25 x 4.9 mm
SMOOTH_VOX = (1.0, 1.0, 1.2)   # (u, v, s) voxels: 0.1, 0.1, 0.6 mm. Along s 0.2 mm left
                               # 0.5 mm ridges where boundaries drift between levels
ISLAND_FRACTION = 0.02
OUTSIDE = 36                   # label index of "not cord"
END_EXT_MM = 3.0               # grid runs this far past the centreline's ends


def write_stl(path, verts, faces):
    tri = verts[faces]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    n /= np.maximum(np.linalg.norm(n, axis=1), 1e-12)[:, None]
    rec = np.zeros(len(faces), dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
    rec["n"], rec["v"] = n, tri
    with open(path, "wb") as f:
        f.write(b"tract body".ljust(80, b" "))
        f.write(struct.pack("<I", len(faces)))
        f.write(rec.tobytes())


def end_planes(paths, S, P, T):
    """The cord's two end caps as planes (point, normal pointing INTO the
    cord), fitted to the STL triangles that face out of each end."""
    tris = []
    for p in paths:
        v, t = read_stl(p)
        tris.append(v[t])
    tri = np.concatenate(tris)
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    area = np.linalg.norm(n, axis=1)/2
    n /= np.maximum(2*area, 1e-12)[:, None]
    cen = tri.mean(1)
    out = []
    for endp, out_dir in ((P[0], -T[0]), (P[-1], T[-1])):
        sel = (n @ out_dir > 0.9) & (np.linalg.norm(cen - endp, axis=1) < 25)
        assert sel.sum() > 10, "no end cap found"
        N = (n[sel]*area[sel, None]).sum(0)
        N /= np.linalg.norm(N)
        p0 = (cen[sel]*area[sel, None]).sum(0)/area[sel].sum()
        assert np.abs((cen[sel] - p0) @ N).max() < 0.05, "end cap is not planar"
        out.append((p0, -N))
    return out


def nearest_fill(lab, unknown, valid):
    """Give unknown voxels the label of the nearest valid voxel (mm metric)."""
    _, idx = distance_transform_edt(~valid, sampling=(GRID_MM, GRID_MM, DS_MM), return_indices=True)
    out = lab.copy()
    out[unknown] = lab[tuple(i[unknown] for i in idx)]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas", default=os.path.join(ROOT, "fem", "out", "atlas"))
    ap.add_argument("--geometry", default=os.path.join(ROOT, "fem", "out", "geometry", "rado_up4_conformal"))
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    out = args.out or os.path.join(args.atlas, "tract_bodies")
    os.makedirs(out, exist_ok=True)
    t0 = time.time()

    maps, meta = {}, {}
    for L in LEVELS:
        d = np.load(os.path.join(args.atlas, "tracts_%s" % L, "tract_map.npz"))
        g = d["u"]
        iu = np.round((np.arange(-U_HALF, U_HALF + 1e-9, GRID_MM) - g[0])/(g[1] - g[0])).astype(int)
        iv = np.round((np.arange(-V_HALF, V_HALF + 1e-9, GRID_MM) - g[0])/(g[1] - g[0])).astype(int)
        maps[L] = d["prob"][:, iu][:, :, iv].astype(np.float32)
        rado = (d["white"] | d["grey"])[iu][:, iv]
        with open(os.path.join(args.atlas, "tracts_%s" % L, "tracts.json")) as f:
            meta[L] = json.load(f)
    u = g[iu]; v = g[iv]
    tracts = meta["T10"]["tracts"]
    with open(os.path.join(ROOT, meta["T10"]["centreline"])) as f:
        cl = json.load(f)
    S, P, T = (np.asarray(cl[k], float) for k in ("s_mm", "points_mm", "tangents"))
    s = np.arange(S[0] - END_EXT_MM, S[-1] + END_EXT_MM + DS_MM/2, DS_MM)
    LEVEL_S = scs_levels.stl_levels(args.geometry, S, P, read_stl)
    assert set(LEVEL_S) == set(LEVELS), LEVEL_S
    print("levels (mm along the cord):", LEVEL_S)

    c = np.column_stack([np.interp(s, S, P[:, d]) for d in range(3)])
    t = np.column_stack([np.interp(s, S, T[:, d]) for d in range(3)])
    lo, hi = s < S[0], s > S[-1]                     # straight on past the ends
    c[lo] = P[0] + (s[lo] - S[0])[:, None]*T[0]
    c[hi] = P[-1] + (s[hi] - S[-1])[:, None]*T[-1]
    t /= np.linalg.norm(t, axis=1)[:, None]
    U = np.array([1.0, 0, 0]) - t[:, :1]*t
    U /= np.linalg.norm(U, axis=1)[:, None]
    V = np.cross(t, U)

    # 1. the anatomy's cord on the grid
    UU, VV = np.meshgrid(u, v, indexing="ij")
    pts = (c[None, None] + UU[:, :, None, None]*U[None, None] + VV[:, :, None, None]*V[None, None])
    pts = pts.reshape(-1, 3)
    stl = lambda pat: glob.glob(os.path.join(args.geometry, pat))
    paths = stl("*whitemater*") + stl("*GreyMater*")
    assert len(paths) == 2, paths
    cord = np.zeros(len(pts), bool)
    for p in paths:
        cord |= InsideTester(*read_stl(p))(pts)
    shape = (len(u), len(v), len(s))
    cord = cord.reshape(shape)
    print("cord voxels %d (%.1f mm3) in %.0f s" % (cord.sum(), cord.sum()*GRID_MM**2*DS_MM, time.time() - t0))
    # carry the full section out past both ends; the end planes trim it (step 6)
    area = cord.sum(axis=(0, 1))
    full = np.flatnonzero(area >= 0.98*np.median(area[area > 0]))
    cord_ext = cord.copy()
    cord_ext[:, :, :full[0]] = cord[:, :, full[0]][:, :, None]
    cord_ext[:, :, full[-1] + 1:] = cord[:, :, full[-1]][:, :, None]
    planes = end_planes(paths, S, P, T)
    print("end planes:", [(p0.round(2).tolist(), nn.round(4).tolist()) for p0, nn in planes])

    # 2. most probable tract, blended between vertebrae
    ks = np.array([LEVEL_S[L] for L in LEVELS])
    lab = np.full(shape, OUTSIDE, np.int16)
    known = np.zeros(shape, bool)
    for j, sj in enumerate(s):
        k = np.clip(np.searchsorted(ks, sj) - 1, 0, len(ks) - 2)
        a = np.clip((sj - ks[k])/(ks[k + 1] - ks[k]), 0, 1)
        p = (1 - a)*maps[LEVELS[k]] + a*maps[LEVELS[k + 1]]
        lab[:, :, j] = np.where(rado, p.argmax(0), OUTSIDE)
        known[:, :, j] = rado
    lab = np.where(cord_ext, lab, OUTSIDE)
    valid = cord_ext & known
    lab = nearest_fill(lab, cord_ext & ~known, valid)

    # 3. islands
    moved = 0
    for i in range(36):
        comp, n = cc_label(lab == i)
        if n <= 1:
            continue
        size = np.bincount(comp.ravel())[1:]
        small = np.isin(comp, 1 + np.flatnonzero(size < ISLAND_FRACTION*size.max()))
        if small.any():
            moved += int(small.sum())
            ok = cord_ext & ~small & (lab != i)
            lab = nearest_fill(lab, small, ok)
    print("islands relabelled: %d voxels" % moved)

    # 4. shared surfaces from smoothed indicators
    present = [i for i in range(36) if (lab == i).any()]
    labels = present + [OUTSIDE]
    pad = lambda a, val: np.pad(a, 1, constant_values=val)
    sm = np.stack([gaussian_filter(pad((lab == k).astype(np.float32), 1.0 if k == OUTSIDE else 0.0), SMOOTH_VOX)
                   for k in labels])
    order = np.argsort(sm, axis=0)
    top1, top2 = order[-1], order[-2]
    m1 = np.take_along_axis(sm, top1[None], 0)[0]
    m2 = np.take_along_axis(sm, top2[None], 0)[0]

    index = []
    for n_, i in enumerate(present):
        tr = tracts[i]
        f = sm[n_] - np.where(top1 == n_, m2, m1)
        if f.max() <= 0:
            print("%-55s skipped: vanishes after smoothing" % tr["name"])
            continue
        verts, faces, _, _ = marching_cubes(f, 0.0)
        verts -= 1                                   # undo the pad
        iu_, iv_, js = verts.T
        uu, vv = u[0] + iu_*GRID_MM, v[0] + iv_*GRID_MM
        j0 = np.clip(np.floor(js).astype(int), 0, len(s) - 2)
        w = np.clip(js - j0, 0, 1)[:, None]
        lerp = lambda A: (1 - w)*A[j0] + w*A[j0 + 1]
        xyz = lerp(c) + uu[:, None]*lerp(U) + vv[:, None]*lerp(V)
        name = "%02d_%s_%s.stl" % (i, tr["side"], tr["short"].split(" (")[0].replace(" ", "_"))
        tri = np.array(faces[:, ::-1], np.uint32, order="C")            # outward
        # writable C-contiguous copies: manifold3d's binding rejects read-only views
        body = m3d.Manifold(m3d.Mesh(vert_properties=np.array(xyz, np.float32, order="C"), tri_verts=tri))
        trimmed = body.status() == m3d.Error.NoError
        if trimmed:
            for p0, nn in planes:
                body = body.trim_by_plane(tuple(nn), float(nn @ p0))
            if body.is_empty():
                print("%-55s skipped: nothing inside the cord's ends" % tr["name"])
                continue
            mesh = body.to_mesh()
            xyz, tri = np.asarray(mesh.vert_properties)[:, :3], np.asarray(mesh.tri_verts)
        else:
            print("%-55s NOT trimmed: %s" % (name, body.status()))
        write_stl(os.path.join(out, name), xyz, tri)
        in_cord = (lab == i) & cord
        nvox = int(in_cord.sum())
        index.append({**{k: tr[k] for k in ("id", "name", "short", "side", "tissue", "colour")},
                      "stl": name, "triangles": int(len(tri)), "ends_trimmed": bool(trimmed),
                      "volume_mesh_mm3": round(float(body.volume()), 2) if trimmed else None,
                      "volume_vox_mm3": round(nvox*GRID_MM**2*DS_MM, 2),
                      "mean_area_mm2": round(nvox*GRID_MM**2/max(1, in_cord.any(axis=(0, 1)).sum()), 3)})
        print("%-55s %7d tris  %.2f mm2" % (name, len(tri), index[-1]["mean_area_mm2"]))
    for old in glob.glob(os.path.join(out, "*.stl")):
        if os.path.basename(old) not in {b["stl"] for b in index}:
            os.remove(old)
    cs = np.flatnonzero(cord.any(axis=(0, 1)))
    with open(os.path.join(out, "bodies.json"), "w") as f:
        json.dump({"levels_s_mm": LEVEL_S, "levels_definition": "vertebra volume centroid projected on the centreline (scs_levels.py)", "ds_mm": DS_MM, "grid_mm": GRID_MM, "smooth_vox": SMOOTH_VOX,
                   "island_fraction": ISLAND_FRACTION, "islands_relabelled_vox": moved,
                   "cord_s_range_mm": [float(s[cs[0]]), float(s[cs[-1]])],
                   "cord_volume_vox_mm3": round(float(cord.sum())*GRID_MM**2*DS_MM, 1),
                   "end_planes": [{"point_mm": p0.round(4).tolist(), "inward_normal": nn.round(6).tolist()}
                                  for p0, nn in planes],
                   "bodies_volume_mesh_mm3": round(sum(b["volume_mesh_mm3"] or 0 for b in index), 1),
                   "geometry": [os.path.relpath(p, ROOT) for p in paths],
                   "frame": "FIBRE_CENTRELINE frame, as register_tracts.py",
                   "bodies": index, "method": __doc__}, f, indent=1)
    print("done in %.0f s" % (time.time() - t0))


if __name__ == "__main__":
    main()
