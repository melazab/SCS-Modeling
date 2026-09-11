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
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
from stlio import read_stl

# Size-field control, all in mm.
H_MIN = 0.20
H_MAX = 2.50
FIELDS = [
    # (name, h_near, plateau_mm, growth)  ->  h = h_near + growth*max(0, d-plateau)
    ("dura",  0.25, 0.40, 1.0),
    ("lead",  0.30, 0.80, 0.8),
    ("white", 0.50, 1.00, 0.6),
    ("canal", 0.70, 1.00, 0.5),
]
GRID = 0.25       # size-field sample spacing (cached in fem/out/sizefield.npz)
MARGIN = 1.0      # box margin around the epidural envelope


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


def main():
    os.makedirs(C.OUT, exist_ok=True)
    t0 = time.time()

    meshes = {k: read_stl(p) for k, p in C.TISSUE_STL.items()}
    lead_pts = []
    for p in list(C.CONTACT_STL.values()) + [C.INSULATOR_STL]:
        v, t = read_stl(p)
        lead_pts.append(densify(v, t, 0.25))
    clouds = {
        "dura":  densify(*meshes["dura"], 0.45),
        "white": densify(*meshes["white"], 0.70),
        "canal": densify(*meshes["epidural"], 0.90),
        "lead":  np.concatenate(lead_pts),
    }
    print("size-field clouds:", {k: len(v) for k, v in clouds.items()}, flush=True)

    lo = meshes["epidural"][0].min(0) - MARGIN
    hi = meshes["epidural"][0].max(0) + MARGIN
    n = np.maximum(np.ceil((hi - lo) / GRID).astype(int) + 1, 2)
    axes = [lo[i] + GRID * np.arange(n[i]) for i in range(3)]
    print("size grid", n, "=", int(np.prod(n)), "pts", flush=True)

    cache = os.path.join(C.OUT, "sizefield.npz")
    if os.path.exists(cache):
        cd = np.load(cache)
        if np.allclose(cd["lo"], lo) and np.array_equal(cd["n"], n):
            h = cd["h"]
            print("size field reused from cache", flush=True)
        else:
            h = None
    else:
        h = None
    if h is None:
        h = np.full(tuple(n), H_MAX, dtype=np.float32)
        gx, gy, gz = np.meshgrid(*axes, indexing="ij")
        P = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], 1)
        for name, hnear, plateau, growth in FIELDS:
            tree = cKDTree(clouds[name])
            d = np.empty(len(P), dtype=np.float64)
            for s_ in range(0, len(P), 1_000_000):
                d[s_:s_ + 1_000_000] = tree.query(P[s_:s_ + 1_000_000], workers=-1)[0]
            cand = (hnear + growth * np.maximum(0.0, d - plateau)).astype(np.float32)
            h = np.minimum(h, cand.reshape(n))
            print("  field %-6s done %.0fs" % (name, time.time() - t0), flush=True)
        del P, gx, gy, gz
        np.clip(h, H_MIN, H_MAX, out=h)
        np.savez_compressed(cache, h=h, lo=lo, n=n, grid=GRID)
    np.clip(h, H_MIN, H_MAX, out=h)

    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("scs_stripped")
    gmsh.model.occ.addBox(*lo, *(hi - lo))
    gmsh.model.occ.synchronize()

    sf = os.path.join(C.OUT, "sizefield.dat")
    with open(sf, "w") as fh:
        fh.write("%.9g %.9g %.9g\n" % tuple(lo))
        fh.write("%.9g %.9g %.9g\n" % (GRID, GRID, GRID))
        fh.write("%d %d %d\n" % tuple(int(v) for v in n))
        h.astype(np.float64).ravel(order="C").tofile(fh, sep="\n", format="%.6g")
        fh.write("\n")
    fid = gmsh.model.mesh.field.add("Structured")
    gmsh.model.mesh.field.setString(fid, "FileName", sf)
    gmsh.model.mesh.field.setNumber(fid, "TextFormat", 1)
    gmsh.model.mesh.field.setNumber(fid, "SetOutsideValue", 1)
    gmsh.model.mesh.field.setNumber(fid, "OutsideValue", H_MAX)
    gmsh.model.mesh.field.setAsBackgroundMesh(fid)

    gmsh.option.setNumber("Mesh.Algorithm3D", 10)      # HXT: fast parallel Delaunay
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.model.mesh.generate(3)
    print("meshed %.0fs" % (time.time() - t0), flush=True)

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(os.path.join(C.OUT, "mesh.msh"))

    ntag, ncoord, _ = gmsh.model.mesh.getNodes()
    order = np.argsort(ntag)
    nodes = ncoord.reshape(-1, 3)[order]
    remap = np.zeros(int(ntag.max()) + 1, dtype=np.int64)
    remap[ntag[order]] = np.arange(len(ntag))
    etypes, etags, enodes = gmsh.model.mesh.getElements(3)
    tets = remap[np.concatenate([e for ty, e in zip(etypes, enodes) if ty == 4])].reshape(-1, 4)
    gmsh.finalize()

    np.savez_compressed(os.path.join(C.OUT, "mesh.npz"), nodes=nodes, tets=tets,
                        box_lo=lo, box_hi=hi)
    print("nodes %d  tets %d   total %.0fs" % (len(nodes), len(tets), time.time() - t0))


if __name__ == "__main__":
    main()
