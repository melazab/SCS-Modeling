"""Fibre trajectories per tract, for the NEURON threshold runs.

    src/atlas/.venv/bin/python src/atlas/place_fibers.py [--per-side 20] [--out DIR]

Pilot (Mohamed, 2026-10-22): the pain pathways, dorsal columns (PAM50
fasciculus gracilis + cuneatus) as Abeta and the spinothalamic tracts (PAM50
"spinal lemniscus", which lumps the anterior and lateral spinothalamic and the
spinoreticular tracts; the atlas does not separate them) as Adelta.

Each fibre keeps its share of its tract along the cord. The tract maps (per
vertebra, blended along the cord as build_tract_bodies.py does) give the
tract's region every 0.5 mm, kept to its main piece (main_piece). At the reference level the fibres are spread
evenly over it (a centroidal Voronoi layout: each fibre at the centre of the
part of the tract nearest to it); stepping along the cord in both directions,
the layout is updated from the previous step's positions with a few Lloyd
iterations, so each fibre moves only as far as the tract does, and no two
merge. A fibre is placed at the tract pixel nearest its cell's centroid, so it
is always inside the tract, even where the tract is a curved band. Twice as
many candidates as fibres are laid out this way and smoothed over SMOOTH_MM;
those that leave white matter (tested against the model's own white-matter
STL at every step, since the tract maps were sampled on RADO's section at one
level), kink (sideways step over MAX_STEP_MM per 0.5 mm) or spend under
MIN_IN_TRACT of their length in their tract are dropped, and of the rest the
per_side that spread most evenly over the tract at T9 are kept. Polylines are
cut to the cord's own end planes (bodies.json), so every point is inside the
cord.

Tried first, and why not: one fixed position within the part of the tract
present at every level (the spinothalamic tracts move between T6 and T12, so
that part is 0.1-0.2 mm2 of ~2.5); and moving each fibre with the tract's
centroid and shape (points left the tract and, nudged back, piled up on its
edge, two fibres per point).

Writes <out>/<tract>_<side>_<nn>.json in the format src/neuron/run_axon.py
reads (units mm, points in the FreeCAD document frame, allowed_tissues), and
<out>/fibers.json indexing them. The threshold runner takes a window of each
centred on the lead it tests.
"""
import argparse
import json
import os
import sys

import numpy as np
from scipy.ndimage import binary_erosion, binary_fill_holes, gaussian_filter1d, label as cc_label
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
ATLAS = os.path.join(ROOT, "fem", "out", "atlas")
# tract name -> (PAM50 ids per side, axon kind, description)
TRACTS = {
    "dc": ({"left": (0, 2), "right": (1, 3)}, "abeta", "dorsal columns (PAM50 fasciculus gracilis + cuneatus)"),
    "st": ({"left": (12,), "right": (13,)}, "adelta",
           "spinothalamic (PAM50 spinal lemniscus: anterior + lateral spinothalamic and spinoreticular)"),
}
DS_MM = 0.5
CORE_MARGIN_MM = 0.05        # inset from the tract edge (incl. its white-matter edge); 0.1
                             # left the spinothalamic band a pixel or two wide near T12
S0_LEVEL = "T9"              # seeds are spread over the tract here (mid T6-T12)
SMOOTH_MM = 4.0              # fibre paths along the cord
CANDIDATES = 2               # laid out per fibre kept; the filter below picks
MAX_STEP_MM = 0.1            # sideways per DS_MM (~11 deg): no kinks
MIN_IN_TRACT = 0.9           # of a fibre's length
END_MARGIN_MM = 0.2


def main_piece(mask, keep=0.1):
    """A tract's region without crumbs: connected pieces under keep x the
    largest are dropped, holes filled. The raw per-voxel winner breaks into
    islands where tracts meet (the spinal lemniscus into up to 9 near T12) and
    fibres hopped onto them; keeping only the largest piece failed the other
    way, where the band, clipped to white matter, splits into two similar
    pieces and "largest" switches between them along the cord."""
    lab, n = cc_label(mask)
    if n > 1:
        size = np.bincount(lab.ravel())[1:]
        mask = np.isin(lab, 1 + np.flatnonzero(size >= keep*size.max()))
    return binary_fill_holes(mask)


def farthest_points(P, n, seed=0):
    rng = np.random.default_rng(seed)
    pick = [int(rng.integers(len(P)))]
    d = np.linalg.norm(P - P[pick[0]], axis=1)
    for _ in range(n - 1):
        pick.append(int(d.argmax()))
        d = np.minimum(d, np.linalg.norm(P - P[pick[-1]], axis=1))
    return P[pick]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-side", type=int, default=20)
    ap.add_argument("--out", default=os.path.join(ATLAS, "fibers", "pilot"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    with open(os.path.join(ATLAS, "tract_bodies", "bodies.json")) as f:
        bodies = json.load(f)
    levels = bodies["levels_s_mm"]
    order = sorted(levels, key=lambda L: levels[L])          # caudal (T12) first
    ks = np.array([levels[L] for L in order])
    maps = {}
    for L in order:
        d = np.load(os.path.join(ATLAS, "tracts_%s" % L, "tract_map.npz"))
        maps[L] = d["prob"].astype(np.float32)
        g, white = d["u"], d["white"]
    h = g[1] - g[0]
    with open(os.path.join(ATLAS, "tracts_T10", "tracts.json")) as f:
        cl_path = json.load(f)["centreline"]
    with open(os.path.join(ROOT, cl_path)) as f:
        cl = json.load(f)
    S, P, T = (np.asarray(cl[k], float) for k in ("s_mm", "points_mm", "tangents"))
    planes = [(np.array(p["point_mm"]), np.array(p["inward_normal"])) for p in bodies["end_planes"]]
    s = np.arange(np.ceil(bodies["cord_s_range_mm"][0]), bodies["cord_s_range_mm"][1] + 1e-9, DS_MM)

    # labels along the cord (blended between vertebrae, constant past T6/T12)
    labels = []
    for sj in s:
        k = int(np.clip(np.searchsorted(ks, sj) - 1, 0, len(ks) - 2))
        a = float(np.clip((sj - ks[k])/(ks[k + 1] - ks[k]), 0, 1))
        p = (1 - a)*maps[order[k]] + a*maps[order[k + 1]]
        labels.append(np.where(white, p.argmax(0), -1).astype(np.int8))
    labels = np.array(labels)

    c = np.column_stack([np.interp(s, S, P[:, d]) for d in range(3)])
    t = np.column_stack([np.interp(s, S, T[:, d]) for d in range(3)])
    t /= np.linalg.norm(t, axis=1)[:, None]
    U = np.array([1.0, 0, 0]) - t[:, :1]*t
    U /= np.linalg.norm(U, axis=1)[:, None]
    V = np.cross(t, U)
    UU, VV = np.meshgrid(g, g, indexing="ij")
    erode = int(round(CORE_MARGIN_MM/h))

    # The tract maps were sampled on RADO's own section at one level; the
    # extended model's white matter sits a little differently along the cord,
    # and the spinothalamic tracts border both the grey matter and the cord's
    # surface (4 % of their points fell outside white matter). So candidate
    # pixels are tested against the model's white-matter STL at every step.
    sys.path.insert(0, os.path.join(ROOT, "fem", "scripts"))
    from inside import InsideTester
    from stlio import read_stl
    wm_stl = [os.path.join(ROOT, p) for p in bodies["geometry"] if "whitemater" in p][0]
    in_wm = InsideTester(*read_stl(wm_stl))
    cand = np.zeros(labels.shape[1:], bool)
    for ids_by_side, _k, _d in TRACTS.values():
        for ids in ids_by_side.values():
            cand |= np.isin(labels, ids).any(axis=0)
    ci = np.nonzero(cand)
    white3d = np.zeros(labels.shape, bool)
    for j in range(len(s)):
        q = c[j] + UU[ci][:, None]*U[j] + VV[ci][:, None]*V[j]
        white3d[j][ci] = in_wm(q)
    print("white matter tested at %d points" % (len(ci[0])*len(s)))

    index, counts = [], {}
    j0 = int(np.argmin(abs(s - levels[S0_LEVEL])))
    for tract, (ids_by_side, kind, desc) in TRACTS.items():
        for side, ids in ids_by_side.items():
            masks = np.array([main_piece(mk) for mk in np.isin(labels, ids) & white3d])
            region = lambda j: np.column_stack([UU[inner[j]], VV[inner[j]]])
            inner = [binary_erosion(mk, iterations=erode) if erode else mk for mk in masks]   # 0 = erode to nothing

            def lloyd(X, centres, iters):
                """Centroidal Voronoi update; each centre ends on the region pixel
                nearest its cell's centroid, so it is inside the region."""
                for _ in range(iters):
                    _, own = cKDTree(centres).query(X)
                    new = centres.copy()
                    for k in range(len(centres)):
                        cell = X[own == k]
                        if len(cell):
                            cen = cell.mean(0)
                            new[k] = cell[np.argmin(((cell - cen)**2).sum(1))]
                    centres = new
                return centres

            npx = np.array([m_.sum() for m_ in inner])
            # a step is "full" when all of the tract lies inside both end planes
            # (a test on its area instead froze the fibres wherever the tract is
            # simply smaller, e.g. the spinothalamic below T12)
            full = np.zeros(len(s), bool)
            for j in range(len(s)):
                if npx[j]:
                    q = c[j] + UU[inner[j]][:, None]*U[j] + VV[inner[j]][:, None]*V[j]
                    full[j] = all(((q - p0) @ nn >= END_MARGIN_MM).all() for p0, nn in planes)
            M = CANDIDATES*args.per_side                     # oversample, then filter
            paths = np.zeros((len(s), M, 2))
            X0 = region(j0)
            paths[j0] = lloyd(X0, farthest_points(X0, M), 30)
            for rng_ in (range(j0 + 1, len(s)), range(j0 - 1, -1, -1)):
                prev = paths[j0]
                for j in rng_:
                    # where an oblique end plane cuts the tract, hold (end_planes cut the fibre)
                    paths[j] = prev = lloyd(region(j), prev, 3) if full[j] else prev
            # Smooth along the cord: where the even layout reorganises (the
            # thin spinothalamic band, T11-T12) a fibre would otherwise step
            # sideways up to 0.6 mm in 0.5 mm, a kink that changes its
            # threshold. Nothing is snapped back (snapping made new kinks and
            # put two fibres on one pixel); fibres that leave white matter, kink
            # or leave their tract are dropped below instead.
            paths = gaussian_filter1d(paths, SMOOTH_MM/DS_MM, axis=0, mode="nearest")
            kept = np.zeros((len(s), M), bool)
            for n in range(M):
                xyz = c + paths[:, n, :1]*U + paths[:, n, 1:]*V
                kept[:, n] = np.all([(xyz - p0) @ nn >= END_MARGIN_MM for p0, nn in planes], axis=0)
            iu = np.clip(np.round((paths[..., 0] - g[0])/h).astype(int), 0, len(g) - 1)
            iv = np.clip(np.round((paths[..., 1] - g[0])/h).astype(int), 0, len(g) - 1)
            in_tract = masks[np.arange(len(s))[:, None], iu, iv]
            in_white = white3d[np.arange(len(s))[:, None], iu, iv]
            step = np.linalg.norm(np.diff(paths, axis=0), axis=2)*(kept[1:] & kept[:-1])
            frac_t = np.array([in_tract[kept[:, n], n].mean() for n in range(M)])
            ok = ((np.array([in_white[kept[:, n], n].all() for n in range(M)]))
                  & (step.max(0) <= MAX_STEP_MM) & (frac_t >= MIN_IN_TRACT))
            good = np.flatnonzero(ok)
            # of those, the per_side that spread most evenly over the tract at T9
            pick = good[np.argsort(-frac_t[good], kind="stable")]
            if len(pick) > args.per_side:
                sel = [int(pick[0])]
                d = np.linalg.norm(paths[j0, pick] - paths[j0, pick[0]], axis=1)
                for _ in range(args.per_side - 1):
                    sel.append(int(pick[d.argmax()]))
                    d = np.minimum(d, np.linalg.norm(paths[j0, pick] - paths[j0, sel[-1]], axis=1))
                pick = np.array(sel)
            pick = pick[np.lexsort((paths[j0, pick, 0], paths[j0, pick, 1]))]      # stable numbering
            for n_out, n in enumerate(pick):
                uu, vv = paths[:, n, 0], paths[:, n, 1]
                xyz = c + uu[:, None]*U + vv[:, None]*V
                keep = kept[:, n]
                name = "%s_%s_%02d" % (tract, side, n_out)
                fib = {"name": name, "units": "mm", "points": np.round(xyz[keep], 4).tolist(),
                       "allowed_tissues": ["white"], "anatomically_validated": False,
                       "purpose": "Pilot fibre: %s, patient %s; keeps its share of the tract along "
                                  "T6-T12 (src/atlas/place_fibers.py)" % (desc, side),
                       "tract": tract, "side": side, "kind": kind,
                       "uv_mm_at_%s" % S0_LEVEL: [round(float(x), 4) for x in paths[j0, n]],
                       "s_range_mm": [float(s[keep][0]), float(s[keep][-1])],
                       "fraction_in_tract": round(float(frac_t[n]), 4),
                       "max_sideways_step_mm": round(float(step[:, n].max()), 4)}
                with open(os.path.join(args.out, name + ".json"), "w") as f:
                    json.dump(fib, f)
                index.append({k: fib[k] for k in ("name", "tract", "side", "kind", "s_range_mm",
                                                  "fraction_in_tract", "max_sideways_step_mm")})
            P = paths[:, pick]
            D = np.linalg.norm(P[:, :, None] - P[:, None], axis=-1)
            D[:, np.arange(len(pick)), np.arange(len(pick))] = np.inf
            D[~(kept[:, pick][:, :, None] & kept[:, pick][:, None, :])] = np.inf
            counts[(tract, side)] = (M, len(good), len(pick))
            print("%-3s %-5s %-6s %d candidates, %d qualify, %d kept; in tract min %.1f %%; spacing "
                  "median %.2f mm, worst %.2f mm; largest sideways step %.3f mm per %.1f mm"
                  % (tract, side, kind, M, len(good), len(pick), 100*frac_t[pick].min(),
                     np.median(D.min(2)[np.isfinite(D.min(2))]), D.min(), step[:, pick].max(), DS_MM))
    with open(os.path.join(args.out, "fibers.json"), "w") as f:
        json.dump({"tracts": {k: {"pam50_ids": v[0], "kind": v[1], "description": v[2]} for k, v in TRACTS.items()},
                   "per_side": args.per_side, "core_margin_mm": CORE_MARGIN_MM, "ds_mm": DS_MM,
                   "seed_level": S0_LEVEL, "smooth_mm": SMOOTH_MM, "candidates_per_fibre": CANDIDATES,
                   "max_sideways_step_mm": MAX_STEP_MM, "min_in_tract": MIN_IN_TRACT,
                   "candidates_qualified_kept": {"%s_%s" % k: v for k, v in counts.items()},
                   "frame": "FIBRE_CENTRELINE frame, as register_tracts.py", "levels_s_mm": levels,
                   "fibers": index, "method": __doc__}, f, indent=1)


if __name__ == "__main__":
    main()
