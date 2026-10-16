"""Make an anatomical STL set fused, smoothly blended and overlap-free.

    src/extend/.venv/bin/python src/extend/conform.py                 # rado_up4 -> rado_up4_conformal
    src/extend/.venv/bin/python src/extend/conform.py --only sympathetic --workers 14

Run it under a memory cap (see meshkit.py). Three steps, per tissue in the FEM's
priority order (fem/scripts/config.py ANATOMY_CLASS, first wins):

1. FUSE. Every body of the tissue is unioned exactly (manifold3d). Pieces that
   touch become one closed body: a sympathetic chain per side, a root complex
   per root, the dural sac with its root sleeves. The sympathetic trunk is first
   rebuilt as one smooth tube (chain_demo.trunk_path), because RADO's straight
   segments meet at angles and offsets that a union would keep as steps.
2. BLEND. A union keeps a sharp crease wherever two pieces meet. The creases
   are found exactly (vertices whose faces came from two different input
   bodies), and only the 4 mm tiles of a global lattice that lie within
   reach of a crease are re-surfaced: signed distance per piece, polynomial
   smooth minimum of width k, marching cubes, clipped to the tile, spliced
   back into the union. Tiles share the global lattice, so neighbouring tiles
   evaluate identical points on their common faces and meet without seams;
   they run in parallel. Away from the joints the union is untouched.
3. PARTITION. Each tissue keeps only the space no higher-priority tissue has
   claimed. Where two tissues overlapped, the lower one is cut back to the
   higher one's surface, so they meet on one shared surface. Because the
   priority is the FEM's own, the result labels space as fem/ already labels
   the original set, except inside the fillets.

Not changed: gaps that exist in RADO (tissues that never touched) stay gaps.
Interfaces coincide geometrically; each side keeps its own triangulation.
"""
import argparse
import collections
import json
import os
import time

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

import meshkit as K
from chain_demo import Grid, exact_distance, point_triangle_distance, smin, surface_samples, trunk_path  # noqa: F401


SIMPLIFY_MM = 0.005
CUT_SIMPLIFY_MM = 1e-5
MIN_PIECE_MM3 = 1e-3
# Loose pieces thinner than this on average are cutting slivers. Measured on
# rado_up4 (2026-10-15), every such piece is <= 40 um; the thinnest real ones
# are 98 um (a trunk fragment), then 140-160 um vessels and >= 430 um roots.
MIN_THICKNESS_MM = 0.05
END_SLAB_MM = 0.3                 # exact (unblended) geometry this deep at the canal's cut ends
END_EXTEND_MM = 0.01              # canal tissues are cut this far past their end plane, then trimmed back


def log(msg):
    print(msg, flush=True)


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------
class Body:
    """An input body: its file, triangles and manifold with a unique original id."""

    def __init__(self, name, tri, M):
        self.name, self.tri, self.M = name, tri, M.as_original()
        self.id = self.M.original_id()
        self.bb = np.array(self.M.bounding_box())
        self._tester = None

    def tester(self):
        if self._tester is None:
            self._tester = K.InsideTester(*K.stlio.weld(self.tri.astype(np.float64)))
        return self._tester


def load(src, only=None):
    by = collections.defaultdict(list)
    for fn in sorted(os.listdir(src)):
        if not fn.lower().endswith(".stl"):
            continue
        pipe = K.pipe_of(os.path.splitext(fn)[0])
        if not pipe or (only and pipe not in only):
            continue
        if pipe == "sympathetic" and "sympathetic_a" in fn:
            continue                      # trunk segments: replaced by one smooth tube per side
        M, _ = K.load_manifold(os.path.join(src, fn))
        by[pipe].append(Body(fn, K.read_tri(os.path.join(src, fn)), M))
    if "sympathetic" in by:
        import chain_demo
        chain_demo.SRC = src
        for side, files in chain_demo.chain_files().items():
            path, r, _, info = trunk_path(files)
            V, T = K.tube_mesh(path, r)
            name = "T8-10 - neuro_sympathetic_a-trunk_%s.STL" % side
            by["sympathetic"].append(Body(name, V[T], K.manifold(V, T)))
            log("  sympathetic trunk %s rebuilt as one tube: %s" % (side, info))
    return by


# --------------------------------------------------------------------------
# Fuse
# --------------------------------------------------------------------------
def boxes_overlap(a, b, pad=0.0):
    return bool(np.all(a[:3] - pad <= b[3:]) and np.all(b[:3] - pad <= a[3:]))


def union_clusters(bodies):
    """Exact union per cluster of bodies that actually intersect. Clustering by
    bounding boxes alone chained all 263 nerve-root bodies into one 89 x 27 x
    227 mm union (adjacent levels' boxes overlap), which ran out of memory;
    by real intersection they are 32 root complexes, fused in 38 s at 0.8 GB.
    Zero-volume slivers are cleaned once per cluster; real voids are kept (they
    may hold another tissue)."""
    import manifold3d as m3
    n = len(bodies)
    rows, cols = [], []
    for i in range(n):
        for j in range(i + 1, n):
            if boxes_overlap(bodies[i].bb, bodies[j].bb, 1e-3) and (bodies[i].M ^ bodies[j].M).volume() > 1e-6:
                rows.append(i)
                cols.append(j)
    _, lab = connected_components(coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n)), directed=False)
    out = []
    for g in np.unique(lab):
        idx = sorted(np.flatnonzero(lab == g), key=lambda i: -bodies[i].M.volume())
        acc = bodies[idx[0]].M
        for i in idx[1:]:
            acc = acc + bodies[i].M
        out.append(K.drop_slivers(acc))
    return m3.Manifold.compose(out)


def crease_points(U, bodies, min_angle_deg=20.0):
    """Points on the creases of a union: edges whose two faces came from
    different input bodies AND meet at an angle of at least min_angle_deg.
    The angle test matters: RADO's layered bodies (a root core in its thin
    middle layer) have nearly coincident surfaces that interpenetrate, so the
    union surface switches body every few mm along a smooth tube. Those
    switches are not creases and must not be blended."""
    mesh = U.to_mesh()
    tv = np.asarray(mesh.tri_verts).astype(np.int64)
    V = np.asarray(mesh.vert_properties)[:, :3].astype(np.float64)
    runs = np.asarray(mesh.run_index)
    oid = np.asarray(mesh.run_original_id)
    face_id = np.repeat(oid, np.diff(runs) // 3)      # run_index counts triangle corners
    known = np.isin(face_id, [b.id for b in bodies])
    n = np.cross(V[tv[:, 1]] - V[tv[:, 0]], V[tv[:, 2]] - V[tv[:, 0]])
    n /= np.maximum(np.linalg.norm(n, axis=1), 1e-30)[:, None]
    e = np.concatenate([tv[:, [0, 1]], tv[:, [1, 2]], tv[:, [2, 0]]])
    f = np.tile(np.arange(len(tv)), 3)
    key = np.sort(e, axis=1)
    order = np.lexsort(key.T[::-1])
    key, f = key[order], f[order]
    same = np.all(key[1:] == key[:-1], axis=1)
    i = np.flatnonzero(same)
    fa, fb = f[i], f[i + 1]
    cross_body = (face_id[fa] != face_id[fb]) & known[fa] & known[fb]
    ang = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", n[fa], n[fb]), -1, 1)))
    sel = cross_body & (ang >= min_angle_deg)
    return V[np.unique(key[i[sel]])]


def crease_tiles(P, margin, tile):
    """Integer indices of the lattice tiles (cubes of side `tile`, on a global
    lattice) within `margin` of any crease point."""
    tiles = set()
    for d in np.array(np.meshgrid(*[[-margin, 0.0, margin]] * 3)).reshape(3, -1).T:
        tiles.update(map(tuple, np.floor((P + d) / tile).astype(int)))
    return sorted(tiles)


# --------------------------------------------------------------------------
# Blend
# --------------------------------------------------------------------------
class Lattice:
    """A box of points of the GLOBAL lattice h * (i, j, k), i, j, k integers.
    Neighbouring tiles therefore evaluate bit-identical points on the planes
    they share, so their re-surfaced patches meet exactly."""

    def __init__(self, i0, i1, h):
        self.i0, self.h = np.asarray(i0), h
        self.shape = tuple(int(x) for x in np.asarray(i1) - self.i0 + 1)
        self.n = int(np.prod(self.shape))
        self.lo = self.i0 * h

    def pts(self):
        ii = np.indices(self.shape).reshape(3, -1).T + self.i0
        return ii * self.h


def body_field(body, pts, band):
    """Signed distance to one body at points, clamped to +-band: exact within
    band of the surface, sign from the body's own inside test."""
    lo, hi = pts.min(0), pts.max(0)
    tri = body.tri
    near = np.all(tri.min(1) <= hi + band, 1) & np.all(tri.max(1) >= lo - band, 1)
    F = np.full(len(pts), band, np.float32)
    if near.any():
        sub = tri[near]
        P, owner = surface_samples(sub, 0.015)
        tree = cKDTree(P)
        u, _ = tree.query(pts, distance_upper_bound=band)
        m = u < band
        if m.any():
            u[m] = exact_distance(pts[m], sub, tree, owner)
        F = np.minimum(u, band).astype(np.float32)
    return np.where(body.tester()(pts), -F, F)


_TILE_JOB = None


def blend_tile(t):
    """Smooth-min union of the class bodies near one lattice tile, re-surfaced
    with marching cubes and clipped to the tile. Returns (verts, tris) or None."""
    import manifold3d as m3
    from skimage.measure import marching_cubes
    bodies, k, h, tile, crease = _TILE_JOB
    step = int(round(tile / h))
    i0 = np.array(t) * step - 3
    lat = Lattice(i0, i0 + step + 6, h)
    pts = lat.pts()
    band = k + 3 * h
    glo, ghi = pts.min(0) - band, pts.max(0) + band
    D, Dmin = None, None
    for b in bodies:                       # fixed order: smin is not associative
        if np.all(b.bb[:3] <= ghi) and np.all(b.bb[3:] >= glo):
            F = body_field(b, pts, band)
            D = F if D is None else smin(D, F, k)
            Dmin = F if Dmin is None else np.minimum(Dmin, F)
    if D is not None:
        # Blend only near creases. The smooth minimum also pushes a surface out
        # by up to k/4 wherever two same-tissue surfaces run within k of each
        # other (a root core inside its thin middle layer: +5% root volume).
        # So: plain union (Dmin) everywhere, the blend within k of a crease
        # point, a smoothstep between k and 2k.
        r1, r2 = k, 2 * k
        cm = np.all((crease >= pts.min(0) - r2) & (crease <= pts.max(0) + r2), 1)
        if cm.any():
            dc, _ = cKDTree(crease[cm]).query(pts, distance_upper_bound=r2)
            u = np.clip((r2 - np.minimum(dc, r2)) / (r2 - r1), 0, 1)
            D = Dmin + (u * u * (3 - 2 * u)).astype(np.float32) * (D - Dmin)
        else:
            D = Dmin
    if D is None:
        return None
    D = D.reshape(lat.shape)
    D[0], D[-1], D[:, 0], D[:, -1], D[:, :, 0], D[:, :, -1] = band, band, band, band, band, band
    if D.min() >= 0:                       # nothing inside away from the border
        return None
    verts, faces, _, _ = marching_cubes(D, level=0.0, spacing=(h, h, h))
    verts += lat.lo
    a, b_, c = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    if np.einsum("ij,ij->i", a, np.cross(b_, c)).sum() < 0:
        faces = faces[:, ::-1]
    P = K.manifold(verts, faces)
    if P.status() != m3.Error.NoError:
        P = K.manifold(*K.weld(verts[faces], 1e-7))
    if P.status() != m3.Error.NoError:
        return ("error", str(P.status()))
    lo = np.array(t) * step * h
    P = P ^ m3.Manifold.cube([step * h] * 3).translate(list(lo))
    if P.is_empty():
        return None
    mesh = P.to_mesh()
    return np.asarray(mesh.vert_properties)[:, :3].copy(), np.asarray(mesh.tri_verts).copy()


def end_regions(src, slab=END_SLAB_MM, margin=1.0, beyond=0.5, reach=8.0):
    """The canal's two cut ends (RADO's bottom, and its top carried up by the
    extension), where the canal tissues all stop in one plane. The epidural
    space, the widest of them, gives each plane and its outward normal: its
    cap there is the largest set of coplanar triangles. The region covers the
    caps of EVERY body ending in that plane within `reach` mm of the epidural
    cap's rim (a vessel end can stick out past the epidural outline), grown by
    `margin` in the plane, from `slab` mm inside the plane to `beyond` mm
    outside it. Returns ([(region, outward normal, plane offset)], the set of
    tissues with a cap in either plane)."""
    import manifold3d as m3
    f = [x for x in os.listdir(src) if "EpiduralSpace" in x and x.lower().endswith(".stl")]
    if len(f) != 1:
        return [], set()

    def faces(tri):
        nn = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        a = 0.5 * np.linalg.norm(nn, axis=1)
        return a, nn / np.maximum(2 * a, 1e-30)[:, None]
    tri = K.read_tri(os.path.join(src, f[0])).astype(np.float64)
    a, nu = faces(tri)
    z = tri.mean(1)[:, 2]
    planes = []
    for near in (z > z.max() - 6.0, z < z.min() + 6.0):
        key, inv = np.unique(np.round(nu[near], 3), axis=0, return_inverse=True)
        best = np.argmax(np.bincount(inv.ravel(), weights=a[near]))
        n = key[best] / np.linalg.norm(key[best])
        P = tri[near][inv.ravel() == best].reshape(-1, 3)
        planes.append((n, float(np.median(P @ n)), P.mean(0), float(np.linalg.norm(P - P.mean(0), axis=1).max())))
    caps, pipes = [[] for _ in planes], set()
    for fn in os.listdir(src):
        pipe = K.pipe_of(os.path.splitext(fn)[0]) if fn.lower().endswith(".stl") else None
        if not pipe:
            continue
        t = K.read_tri(os.path.join(src, fn)).astype(np.float64)
        ta, tn = faces(t)
        for i, (n, d, c, rad) in enumerate(planes):
            on = (np.abs(t @ n - d).max(1) < 0.01) & (np.abs(tn @ n) > 0.999)
            on &= np.linalg.norm(t.mean(1) - c, axis=1) < rad + reach
            caps[i].append(t[on].reshape(-1, 3))
            if ta[on].sum() > 0.01:
                pipes.add(pipe)
    out = []
    for (n, d, c, rad), Ps in zip(planes, caps):
        P = np.vstack(Ps)
        r = P - c
        r = r - np.outer(r @ n, n)
        P = P + margin * r / np.maximum(np.linalg.norm(r, axis=1), 1e-9)[:, None]
        out.append((m3.Manifold.hull_points(np.vstack([P - slab * n, P + beyond * n])), n, d))
    return out, pipes


def past_ends(M, ends, e=END_EXTEND_MM):
    """M carried e mm past each canal end plane, inside that end's region: the
    slice of M within 2e of the plane, moved e outward and unioned on. It
    overlaps M by e, so no faces coincide (a union of exactly touching faces
    leaves them in place)."""
    for R, n, off in ends:
        if not boxes_overlap(np.array(M.bounding_box()), np.array(R.bounding_box())):
            continue
        sl = (M ^ R).trim_by_plane(list(n), off - 2 * e).trim_by_plane(list(-n), -off)
        if not sl.is_empty():
            M = M + sl.translate(list(e * n))
    return M


def fuse_and_blend(pipe, bodies, k, h, tile, workers, ends=()):
    import multiprocessing
    import manifold3d as m3
    global _TILE_JOB
    t0 = time.time()
    U = union_clusters(bodies)
    U_exact = U
    t1 = time.time()
    rep = {"bodies": len(bodies), "union_volume_mm3": round(U.volume(), 2), "union_s": round(t1 - t0, 1)}
    crease = crease_points(U, bodies) if k > 0 else np.zeros((0, 3))
    tiles = crease_tiles(crease, 2 * k + 4 * h, tile) if k > 0 else []
    rep["blend_tiles"] = len(tiles)
    if tiles:
        for b in bodies:                   # build inside testers before forking (shared copy-on-write)
            b.tester()
        _TILE_JOB = (bodies, k, h, tile, crease)
        res, tb = [], time.time()
        with multiprocessing.get_context("fork").Pool(workers) as pool:
            for i, r in enumerate(pool.imap(blend_tile, tiles, chunksize=2)):
                res.append(r)
                if (i + 1) % 50 == 0 or i + 1 == len(tiles):
                    log("    %s: %d / %d tiles blended, %.0f s" % (pipe, i + 1, len(tiles), time.time() - tb))
        _TILE_JOB = None
        for b in bodies:
            b._tester = None
        errors = [r for r in res if isinstance(r, tuple) and isinstance(r[0], str)]
        if errors:
            raise RuntimeError("%s: %d blend tiles failed: %s" % (pipe, len(errors), errors[:3]))
        t2 = time.time()
        step = tile
        cubes = [m3.Manifold.cube([step] * 3).translate(list(np.array(t) * step)) for t in tiles]
        cut = cubes[0]
        for c in cubes[1:]:
            cut = cut + c
        cut.num_tri()
        U = U - cut
        U.num_tri()
        for r in res:
            if r is not None:
                U = U + K.manifold(*r)
                U.num_tri()
        U = K.drop_slivers(U)
        rep.update({"blend_s": round(t2 - t1, 1), "splice_s": round(time.time() - t2, 1),
                    "blended_volume_mm3": round(U.volume(), 2)})
        # At the canal's cut ends the exact union goes back in. Marching cubes
        # on the 0.05 mm lattice cannot keep a tilted flat end flat: it left
        # vessel ends up to 0.044 mm short of the top plane, and the partition
        # then filled that with a film of CSF or epidural (2026-10-15).
        restored = 0.0
        for R, _, _ in ends:
            if (U_exact ^ R).is_empty() and (U ^ R).is_empty():
                continue
            U = (U - R) + (U_exact ^ R)
            U.num_tri()
            restored += (U_exact ^ R).volume()
        if ends:
            U = K.drop_slivers(U)
            rep["end_slab_exact_mm3"] = round(restored, 3)
    # Drops the marching-cubes slivers and pinch edges and thins the 0.05 mm
    # patch triangles ~5x; no surface moves more than SIMPLIFY_MM. Done before
    # the partition, so the cut faces copy the higher-priority surfaces exactly.
    U = U.simplify(SIMPLIFY_MM)
    rep["components"] = len(U.decompose())
    rep["triangles"] = U.num_tri()
    return U, rep


# --------------------------------------------------------------------------
# Partition, naming, output
# --------------------------------------------------------------------------
def drop_thin(M, vmin, tmin):
    """Keep the shells (pieces and voids) with |volume| > vmin and mean
    thickness, |volume| / (area / 2), at least tmin."""
    import manifold3d as m3
    keep = [c for c in M.decompose()
            if abs(c.volume()) > vmin and abs(c.volume()) / max(c.surface_area() / 2, 1e-12) >= tmin]
    return m3.Manifold.compose(keep) if keep else m3.Manifold()


def group_shells(M):
    """Split into bodies: each positive shell with the voids it encloses."""
    import manifold3d as m3
    pos, neg = [], []
    for c in M.decompose():
        v = c.volume()
        if v > 1e-6:
            pos.append(c)
        elif v < -1e-6:
            neg.append(c)
    groups = [[p] for p in pos]
    for c in neg:
        cb = np.array(c.bounding_box())
        host = min((g for g in groups if np.all(np.array(g[0].bounding_box())[:3] <= cb[:3])
                    and np.all(np.array(g[0].bounding_box())[3:] >= cb[3:])),
                   key=lambda g: g[0].volume(), default=None)
        if host is not None:
            host.append(c)
    return [m3.Manifold.compose(g) if len(g) > 1 else g[0] for g in groups]


def name_parts(pipe, parts, bodies):
    by_id = {b.id: b for b in bodies}
    names, used = [], collections.Counter()
    for P in parts:
        mesh = P.to_mesh()
        runs, oid = np.asarray(mesh.run_index), np.asarray(mesh.run_original_id)
        cnt = collections.Counter()
        for i, o in enumerate(oid):
            if o in by_id and runs[i + 1] > runs[i]:
                cnt[o] += int(runs[i + 1] - runs[i]) // 3
        members = [by_id[o].name for o, _ in cnt.most_common()]
        if not members:                    # made entirely of blend patch faces
            c = np.array(P.bounding_box())
            members = [min(bodies, key=lambda b: np.linalg.norm(b.bb.reshape(2, 3).mean(0) - c.reshape(2, 3).mean(0))).name]
        stem = os.path.splitext(members[0])[0]
        used[stem] += 1
        name = "%s__conformal%s.STL" % (stem, "" if used[stem] == 1 else "_%d" % used[stem])
        assert K.pipe_of(os.path.splitext(name)[0]) == pipe, name
        names.append((name, members))
    return names


def save_parts(path, parts, names, rep):
    """Stage-1 checkpoint of one tissue: every part's mesh plus its members."""
    V, T, offs = [], [], [[0, 0]]
    for P in parts:
        mesh = P.to_mesh()
        v, t = np.asarray(mesh.vert_properties)[:, :3], np.asarray(mesh.tri_verts)
        V.append(v)
        T.append(t)
        offs.append([offs[-1][0] + len(v), offs[-1][1] + len(t)])
    np.savez(path + ".npz", verts=np.vstack(V), tris=np.vstack(T), offsets=np.array(offs))
    json.dump({"report": rep, "parts": [{"name": n, "members": m} for n, m in names]}, open(path + ".json", "w"), indent=1)


def load_parts(path):
    z = np.load(path + ".npz")
    meta = json.load(open(path + ".json"))
    o = z["offsets"]
    out = []
    for i, pm in enumerate(meta["parts"]):
        v = z["verts"][o[i][0]:o[i + 1][0]]
        t = z["tris"][o[i][1]:o[i + 1][1]]          # saved per part, already local
        out.append((pm["name"], pm["members"], K.manifold(v, t)))
    return out, meta["report"]


def stage1(a, stage_dir, only):
    """Fuse and blend each tissue, one at a time, checkpointing each to disk."""
    ends, _ = end_regions(a.src)
    log("  canal cut ends: %d (exact geometry kept %.2f mm deep there)" % (len(ends), END_SLAB_MM))
    for pipe in K.ORDER:
        if only and pipe not in only:
            continue
        path = os.path.join(stage_dir, pipe)
        if os.path.exists(path + ".json") and not a.redo:
            log("  %-12s stage 1 checkpoint exists, reused" % pipe)
            continue
        by = load(a.src, {pipe})
        if pipe not in by:
            continue
        k = 0.0 if pipe in set(filter(None, a.no_blend.split(","))) else a.k
        U, rep = fuse_and_blend(pipe, by[pipe], k, a.h, a.tile, a.workers, ends)
        rep["k_mm"] = k
        parts = group_shells(U)
        names = name_parts(pipe, parts, by[pipe])
        save_parts(path, parts, names, rep)
        log("  %-12s %s" % (pipe, rep))
        del by, U, parts


def stage2(a, stage_dir, out):
    """Each part loses whatever any higher-priority part overlapping it claims.
    Part by part, so every boolean stays small."""
    t0 = time.time()
    tissues = {}
    for pipe in K.ORDER:
        path = os.path.join(stage_dir, pipe)
        if os.path.exists(path + ".json"):
            tissues[pipe] = load_parts(path)
    # The canal tissues all end flush in the two cut planes, so cutting one by
    # another there pits coincident end faces against each other; manifold3d
    # resolves that differently at each end, and at the top it left the
    # epidural and CSF with zero-thickness lids over the vessel ends (cream
    # wedges in FreeCAD, 2026-10-15). So those tissues are carried
    # END_EXTEND_MM past the planes, cut there, and trimmed back flush: any
    # debris from the coincident faces lies in the part trimmed off.
    ends, end_pipes = end_regions(a.src)
    beyond = [R.trim_by_plane(list(n), off) for R, n, off in ends]
    log("  canal end planes: %d; cut %.2f mm past them, then trimmed: %s" % (len(ends), END_EXTEND_MM, sorted(end_pipes)))
    outputs, report, higher = [], {"src": a.src, "k_mm": a.k, "h_mm": a.h, "priority": K.ORDER, "tissues": {}}, []
    for pipe in K.ORDER:
        if pipe not in tissues:
            continue
        parts, rep = tissues[pipe]
        at_ends = pipe in end_pipes
        vol_in, vol_out, n_out = 0.0, 0.0, 0
        for name, members, P in parts:
            vol_in += P.volume()
            if at_ends:
                P = past_ends(P, ends)
            bb = np.array(P.bounding_box())
            for hb, H in higher:
                if boxes_overlap(bb, hb):
                    P = P - H
                    P.num_tri()
            if at_ends:
                for B in beyond:
                    P = P - B
                P.num_tri()
            # One file per fused part, however many pieces the cut leaves (a
            # root complex is cut by its DRG and by the cord into ~50 pieces;
            # 1,637 separate root files would swamp FreeCAD). Pieces under
            # MIN_PIECE_MM3, or on average thinner than MIN_THICKNESS_MM, are
            # numerical slivers. RADO's touching tissues have surfaces that
            # nearly coincide and cross back and forth by microns, so cutting
            # one by the other leaves skins mostly 1-3 um thick. As loose sheets
            # up to 16 mm long they passed the volume test: 233 in the
            # vertebrae, as streaks on the canal wall, and 126 of the dura's 128
            # pieces (2026-10-15).
            P = drop_thin(P, MIN_PIECE_MM3, MIN_THICKNESS_MM)
            if P.is_empty():
                continue
            # removes the pinch edges and degenerate triangles that cutting
            # along nearly coincident surfaces leaves (CSF inner vs cord);
            # moves no surface more than 10 nm, so shared faces stay shared
            P = P.simplify(CUT_SIMPLIFY_MM)
            # simplify can pinch sub-micron crumbs off (29 in the dura)
            P = drop_thin(P, MIN_PIECE_MM3, MIN_THICKNESS_MM)
            tri = K.manifold_tri(P)
            K.write_stl(os.path.join(out, name), tri, b"src/extend/conform.py")
            au = K.stlio.audit(*K.stlio.weld(K.read_tri(os.path.join(out, name)), tol=1e-9))
            shells = [c.volume() for c in P.decompose()]
            outputs.append({"file": name, "tissue": pipe, "members": members, "volume_mm3": round(P.volume(), 3),
                            "pieces": sum(v > 0 for v in shells), "voids": sum(v < 0 for v in shells),
                            "triangles": len(tri), "open_edges": au["boundary_edges"],
                            "nonmanifold_edges": au["nonmanifold_edges"]})
            vol_out += P.volume()
            n_out += 1
        # this tissue's FUSED parts (before its own cut) are what lower tissues
        # lose; a canal tissue cuts as carried past the end planes
        for _, _, Pf in parts:
            Hc = past_ends(Pf, ends) if at_ends else Pf
            higher.append((np.array(Hc.bounding_box()), Hc))
        rep.update({"fused_volume_mm3": round(vol_in, 2), "final_volume_mm3": round(vol_out, 2), "bodies_out": n_out})
        report["tissues"][pipe] = rep
        log("  %-12s %.1f mm3 fused -> %.1f mm3 after higher-priority tissues, %d bodies (%.0f s)"
            % (pipe, vol_in, vol_out, n_out, time.time() - t0))
    report["outputs"] = outputs
    bad = [o for o in outputs if o["open_edges"] or o["nonmanifold_edges"]]
    report["not_closed"] = [o["file"] for o in bad]
    report["stage2_s"] = round(time.time() - t0, 1)
    json.dump(report, open(os.path.join(out, "conform_report.json"), "w"), indent=1)
    log("wrote %d bodies to %s; %d not closed" % (len(outputs), out, len(bad)))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--src", default=os.path.join(K.C.OUT, "geometry", "rado_up4"))
    ap.add_argument("--out", default=None)
    ap.add_argument("--k", type=float, default=0.25, help="blend width, mm (0 = plain union)")
    ap.add_argument("--h", type=float, default=0.05, help="blend grid spacing, mm")
    ap.add_argument("--tile", type=float, default=4.0, help="blend tile side, mm (a multiple of --h)")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--only", default=None, help="comma-separated fem tissue names: stage 1 for these only")
    ap.add_argument("--redo", action="store_true", help="recompute stage-1 checkpoints that exist")
    ap.add_argument("--no-blend", default="root",
                    help="comma-separated tissues fused without blending. Default root: RADO's root "
                         "bodies are layered (core in a thin middle layer, rootlets side by side), and a "
                         "blend there adds collars and lumps and +0.7%% (k 0.1) to +3.4%% (k 0.25) root "
                         "volume in the rootlet entry zone. Use '' to blend every tissue.")
    ap.add_argument("--stage", choices=["1", "2", "both"], default="both")
    a = ap.parse_args()
    only = set(a.only.split(",")) if a.only else None
    out = a.out or a.src.rstrip("/") + "_conformal" + ("" if a.k == 0.25 else "_k%g" % a.k)
    stage_dir = os.path.join(out, "stage1")
    os.makedirs(stage_dir, exist_ok=True)
    t0 = time.time()
    if a.stage in ("1", "both"):
        log("stage 1: fuse and blend, %s%s" % (a.src, " (only %s)" % sorted(only) if only else ""))
        stage1(a, stage_dir, only)
    if a.stage in ("2", "both"):
        log("stage 2: partition by priority")
        stage2(a, stage_dir, out)
    log("done in %.0f s" % (time.time() - t0))


if __name__ == "__main__":
    main()
