"""Build a single conformal tetrahedral mesh of the stripped SCS domain.

WHY A SINGLE BOX MESH AND NOT A MULTI-BODY ASSEMBLY
---------------------------------------------------
RADO's compartments tile space, but their STL tessellations do NOT match at the
shared interfaces (measured: the dura and epidural shells share 4 vertices out
of ~1800; CSF and white matter share 237 of 2138).  There is therefore no way
to assemble them into a conformal multi-volume mesh without a boolean/remesh
repair step.  Instead we mesh one box that contains everything, graded by
distance to the structures that matter, and assign a tissue to every tetrahedron
by point classification (assign_and_solve.py).  Interfaces are then resolved to
the local element size rather than followed exactly -- that is the headline
approximation of this route and it is why the mesh is refined hard on the dura.

Output: fem/out/mesh.msh (gmsh 2.2 ASCII) + fem/out/mesh.npz
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
import artifacts
from stlio import read_stl

# Full-anatomy defaults, in mm. These supersede the historical canal-only
# field of record. Changing parameters or geometry invalidates the mesh cache.
DEFAULT_PARAMS = dict(
    H_MIN=0.20,
    H_MAX=2.50,
    FIELDS=[
        # (name, h_near, plateau_mm, growth) -> h = h_near + growth*max(0, d-plateau)
        ("dura",  0.25, 0.40, 1.0),
        ("lead",  0.30, 0.80, 0.8),
        ("white", 0.50, 1.00, 0.6),
        ("canal", 0.70, 1.00, 0.5),
        ("fine", 0.30, 0.20, 0.8),
        ("coarse", 3.0, 0.50, 1.0),
    ],
    GRID=0.25,    # size-field sample spacing (cached as <out_dir>/sizefield.npz)
    MARGIN=1.0,   # box margin around all anatomy and the lead
)


def densify(verts, tris, target):
    """Subdivide triangles until every edge < target; return the point cloud."""
    v, t = verts.copy(), tris.copy()
    pts = [v]
    for _ in range(8):
        a, b, c = v[t[:, 0]], v[t[:, 1]], v[t[:, 2]]
        long = np.maximum.reduce([
            np.linalg.norm(b - a, axis=1),
            np.linalg.norm(c - b, axis=1),
            np.linalg.norm(a - c, axis=1)])
        if long.max() < target:
            break
        sel = long >= target
        a, b, c = a[sel], b[sel], c[sel]
        mid = np.concatenate([(a + b) / 2, (b + c) / 2, (c + a) / 2, (a + b + c) / 3])
        pts.append(mid)
        # 4-way split of the selected triangles, re-welded next round
        n0 = len(v)
        newv = np.concatenate([v, (a + b) / 2, (b + c) / 2, (c + a) / 2])
        k = sel.sum()
        iab = n0 + np.arange(k)
        ibc = iab + k
        ica = ibc + k
        ta, tb, tc = t[sel, 0], t[sel, 1], t[sel, 2]
        newt = np.concatenate([
            t[~sel],
            np.stack([ta, iab, ica], 1), np.stack([iab, tb, ibc], 1),
            np.stack([ica, ibc, tc], 1), np.stack([iab, ibc, ica], 1)])
        v, t = newv, newt
    return np.concatenate(pts)


def configure_gmsh_threads(gmsh, threads):
    """Serialize boundary meshing to avoid Structured-field initialization races."""
    threads = max(1, int(threads))
    gmsh.option.setNumber("General.NumThreads", threads)
    # Gmsh 4.15.2 crashes with heap corruption when parallel edge meshing
    # first accesses our Structured background field. Load/use it serially
    # through boundary meshing; retain the requested parallelism for HXT.
    gmsh.option.setNumber("Mesh.MaxNumThreads1D", 1)
    gmsh.option.setNumber("Mesh.MaxNumThreads2D", 1)
    gmsh.option.setNumber("Mesh.MaxNumThreads3D", threads)
    print("Gmsh thread limits: edges/surfaces 1; HXT volume %d "
          "(parallel stages may use fewer)" % threads, flush=True)


def build(params=None, out_dir=None, contact_stl=None, insulator_stl=None, progress_cb=None):
    """Build the mesh with the given size-field `params` (see DEFAULT_PARAMS
    for the schema) into `out_dir`. Both default to the field-of-record
    values, so `build()` with no arguments is exactly what this script has
    always done. Returns (nodes, tets, lo, hi).

    `contact_stl`/`insulator_stl` override which lead geometry the size
    field is graded against (the "lead" region in `params["FIELDS"]") --
    default to `config.CONTACT_STL`/`config.INSULATOR_STL` (the frozen
    fem/leads/fem_dorsal_T10 lead) if omitted. This is an explicit
    parameter, not a mutation of config.py's module state, so a caller
    building one lead's mesh can never leave a stale override behind for
    the next caller in the same process (see
    fem/scripts/SCS_Mesh_Generator.FCMacro's module docstring). Tissue
    geometry (config.TISSUE_STL) is never overridden -- only the lead is
    swappable, everything else is RADO's fixed anatomy.

    `progress_cb`, if given, is called with an int covering the two slow
    parts of this function on a 0-75 SCALE (mesh_preview.py, the one real
    caller, owns the rest of its own 0-100 stdout scale for the
    classify/dura-check/preview-write stages that happen after build()
    returns -- see that module's docstring): the size-field KDTree
    computation (one call per region in `params["FIELDS"]`, landing on 60
    when the loop finishes, or a single jump straight to 60 if the on-disk
    sizefield cache is reused instead of recomputed) then 75 once gmsh's
    own tetrahedralization finishes. Default None (main() below and every
    existing caller pass nothing) so this is purely additive -- no
    behaviour change for the field-of-record path.
    """
    def _progress(pct):
        if progress_cb is not None:
            progress_cb(pct)
    # Imported here, not at module level: this keeps `import build_mesh`
    # (done just to read DEFAULT_PARAMS, e.g. from
    # fem/scripts/SCS_Mesh_Generator.FCMacro running in FreeCAD's OWN Python,
    # which has no scipy) working even where scipy isn't installed. build()
    # itself is only ever actually called from fem/.venv (this script direct,
    # or mesh_preview.py's subprocess), which does have it.
    from scipy.spatial import cKDTree
    p = dict(DEFAULT_PARAMS) if params is None else params
    out = C.OUT if out_dir is None else out_dir
    contact_stl = C.CONTACT_STL if contact_stl is None else contact_stl
    insulator_stl = C.INSULATOR_STL if insulator_stl is None else insulator_stl
    h_min, h_max, fields, grid, margin = (
        p["H_MIN"], p["H_MAX"], p["FIELDS"], p["GRID"], p["MARGIN"])

    os.makedirs(out, exist_ok=True)
    t0 = time.time()

    from mesh_geometry import geometry, grid_points
    signature = artifacts.mesh_signature(p, contact_stl, insulator_stl)
    mesh_path = os.path.join(out, 'mesh.npz')
    if artifacts.valid(mesh_path, signature):
        with np.load(mesh_path) as cached:
            print('Reusing completed volume mesh; resuming classification/preview', flush=True)
            _progress(75)
            return cached['nodes'], cached['tets'], cached['box_lo'], cached['box_hi']
    if h_min <= 0 or h_max < h_min or grid <= 0 or margin < 0:
        raise ValueError('Require 0 < H_MIN <= H_MAX, GRID > 0 and MARGIN >= 0')
    # Enforce the budget here, not only in the GUI wrapper: run_all.sh and
    # direct callers must get the same protection.
    import mesh_cost
    print('Checking full mesh/preview memory budget', flush=True)
    estimate = mesh_cost.CostModel(contact_stl, insulator_stl, margin=margin).estimate(p)
    severity, description = mesh_cost.describe(estimate)
    print(description, flush=True)
    if severity == 'danger':
        raise RuntimeError('Mesh refused before build: ' + description)
    clouds, bounds = geometry(contact_stl, insulator_stl)
    print('size-field clouds:', {k: len(v) for k, v in clouds.items()}, flush=True)
    lo, hi = bounds.min(0) - margin, bounds.max(0) + margin
    n = np.maximum(np.ceil((hi - lo) / grid).astype(int) + 1, 2)
    axes = [lo[i] + grid * np.arange(n[i]) for i in range(3)]
    print("size grid", n, "=", int(np.prod(n)), "pts", flush=True)

    cache_key = signature
    cache = os.path.join(out, "sizefield.npz")
    if os.path.exists(cache):
        cd = np.load(cache, allow_pickle=True)
        # "params" wasn't recorded before build() took a params argument;
        # treat an old cache with no params key as a plain mismatch (rebuild)
        # rather than crashing on a missing key.
        cached_params = cd["params"].item() if "params" in cd.files else None
        if (np.allclose(cd["lo"], lo) and np.array_equal(cd["n"], n)
                and cached_params == cache_key):
            h = cd["h"]
            print("size field reused from cache", flush=True)
        else:
            h = None
    else:
        h = None
    if h is None:
        h = np.full(tuple(n), h_max, dtype=np.float32)
        count = int(np.prod(n))
        chunk = 200_000
        n_chunks = max(1, (count + chunk - 1) // chunk)
        for i, (name, hnear, plateau, growth) in enumerate(fields):
            tree = cKDTree(clouds.pop(name))
            if hnear <= 0 or plateau < 0 or growth <= 0:
                raise ValueError('Region sizes/growth must be positive and plateau nonnegative')
            radius = plateau + max(0, h_max - hnear) / growth + 1e-9
            for j, start in enumerate(range(0, count, chunk)):
                stop = min(count, start + chunk)
                points = grid_points(lo, n, grid, start, stop)
                distance = tree.query(points, workers=int(os.environ.get('OMP_NUM_THREADS', '1')),
                                      distance_upper_bound=radius)[0]
                candidate = hnear + growth * np.maximum(0, distance - plateau)
                np.minimum(h.ravel()[start:stop], candidate, out=h.ravel()[start:stop])
                print('  field %s chunk %d/%d %.0fs' % (name, j + 1, n_chunks, time.time() - t0), flush=True)
                _progress(round(60 * (i + (j + 1) / n_chunks) / len(fields)))
            del tree
        del clouds
        np.clip(h, h_min, h_max, out=h)
        np.savez_compressed(cache, h=h, lo=lo, n=n, grid=grid, params=cache_key)
    else:
        # Cache hit: the slow part above never ran, so jump the progress
        # bar straight past its budget instead of leaving it looking stuck
        # at 0% for a build that is actually already 60% "done".
        _progress(60)
    np.clip(h, h_min, h_max, out=h)

    import gmsh
    gmsh.initialize()
    # OMP_NUM_THREADS alone does not override Gmsh's default of one thread.
    threads = max(1, int(os.environ.get("OMP_NUM_THREADS", "1")))
    configure_gmsh_threads(gmsh, threads)
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("scs_stripped")
    gmsh.model.occ.addBox(*lo, *(hi - lo))
    gmsh.model.occ.synchronize()

    sf = os.path.join(out, "sizefield.dat")
    with open(sf, "w") as fh:
        fh.write("%.9g %.9g %.9g\n" % tuple(lo))
        fh.write("%.9g %.9g %.9g\n" % (grid, grid, grid))
        fh.write("%d %d %d\n" % tuple(int(v) for v in n))
        h.astype(np.float64).ravel(order="C").tofile(fh, sep="\n", format="%.6g")
        fh.write("\n")
    fid = gmsh.model.mesh.field.add("Structured")
    gmsh.model.mesh.field.setString(fid, "FileName", sf)
    gmsh.model.mesh.field.setNumber(fid, "TextFormat", 1)
    gmsh.model.mesh.field.setNumber(fid, "SetOutsideValue", 1)
    gmsh.model.mesh.field.setNumber(fid, "OutsideValue", h_max)
    gmsh.model.mesh.field.setAsBackgroundMesh(fid)

    gmsh.option.setNumber("Mesh.Algorithm3D", 10)      # HXT: fast parallel Delaunay
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.model.mesh.generate(3)
    print("meshed %.0fs" % (time.time() - t0), flush=True)
    _progress(75)

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(os.path.join(out, "mesh.msh"))

    ntag, ncoord, _ = gmsh.model.mesh.getNodes()
    order = np.argsort(ntag)
    nodes = ncoord.reshape(-1, 3)[order]
    remap = np.zeros(int(ntag.max()) + 1, dtype=np.int64)
    remap[ntag[order]] = np.arange(len(ntag))
    etypes, etags, enodes = gmsh.model.mesh.getElements(3)
    tets = remap[np.concatenate([e for ty, e in zip(etypes, enodes) if ty == 4])].reshape(-1, 4)
    gmsh.finalize()

    if signature != artifacts.mesh_signature(p, contact_stl, insulator_stl):
        raise RuntimeError('Model inputs changed during build; mesh not published')
    artifacts.atomic_npz(mesh_path, nodes=nodes, tets=tets, box_lo=lo, box_hi=hi)
    artifacts.publish(mesh_path, signature)
    print("nodes %d  tets %d   total %.0fs" % (len(nodes), len(tets), time.time() - t0))
    return nodes, tets, lo, hi


def main():
    build()


if __name__ == "__main__":
    main()
