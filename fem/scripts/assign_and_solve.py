"""Assign tissues to tets and solve div(sigma grad V) = 0 with a P1 FEM.

PHYSICS
-------
Steady current conduction, exactly the equation in the brief:
    div(sigma grad V) = 0,   J = -sigma grad V
discretised with linear (P1) tetrahedra.  Values are VOLTS, conductivity is
S/m, current is amperes.  Nothing here is an analogy to another physics and
nothing is relabelled.

BOUNDARY CONDITIONS AS ACTUALLY IMPOSED
---------------------------------------
* Current injection is a volumetric source spread over the nodes of the driven
  contact body, weighted by nodal volume.  Because the contact is metal
  (sigma = 1e4 S/m, >2e5x its surroundings) it is an equipotential body, so
  where inside it the current is injected does not affect the exterior field;
  the metal redistributes it over its own surface.  Verified by reporting the
  potential spread within each contact.
* The six undriven contacts are left as plain high-conductivity bodies with no
  source term.  That is an EXACT floating conductor, not an approximation of
  one: with no source inside, conservation forces net current through the body
  to be zero, and the high conductivity makes it equipotential.  Both are
  measured and reported rather than assumed.
* The outer boundary of the tissue domain is insulating.  This is the natural
  ("do nothing") boundary condition of the weak form, so it needs no code: every
  face with no neighbouring element carries zero normal current.
* The all-Neumann problem is singular up to an additive constant.  The gauge is
  fixed by pinning ONE node -- the mesh node furthest from the lead axis -- to
  0 V.  Because the source terms sum to exactly zero the pinned node draws no
  real current; the reaction there is reported as a fraction of 1 A and is the
  check that the pin is a gauge choice and not a current path.

BASIS FIELDS
------------
Eight solves are done, one per contact: +1 A into contact i, -1 A spread over
the outer boundary of the domain, weighted by boundary nodal area.  Any
zero-net-current contact configuration is then a linear combination,
V = sum_i I_i phi_i, with the boundary return cancelling identically.  The
requested bipolar drive is phi_3 - phi_5 at 1 A.
"""
import os
import sys
import time

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
from stlio import read_stl
from inside import InsideTester

INCLUDE_BACKGROUND = os.environ.get("SCS_BACKGROUND", "0") == "1"


def classify_tets(nodes, tets):
    """Label every tet by the tissue containing its centroid (ORDER: first wins)."""
    cen = nodes[tets].mean(axis=1)
    testers = {}
    for i, p in C.CONTACT_STL.items():
        testers["contact%d" % i] = InsideTester(*read_stl(p))
    testers["insulator"] = InsideTester(*read_stl(C.INSULATOR_STL))
    for k, p in C.TISSUE_STL.items():
        testers[k] = InsideTester(*read_stl(p))
    lab = np.full(len(cen), -1, dtype=np.int16)
    todo = np.arange(len(cen))
    for i, name in enumerate(C.ORDER):
        if todo.size == 0:
            break
        hit = np.zeros(todo.size, dtype=bool)
        for s in range(0, todo.size, 400_000):
            sl = slice(s, s + 400_000)
            hit[sl] = testers[name](cen[todo[sl]])
        lab[todo[hit]] = i
        todo = todo[~hit]
    return lab


def sigma_of(lab):
    s = np.full(len(lab), C.SIGMA["background"] if INCLUDE_BACKGROUND else np.nan)
    for i, name in enumerate(C.ORDER):
        v = C.SIGMA_METAL if name.startswith("contact") else C.SIGMA[name]
        s[lab == i] = v
    return s


def assemble(nodes, tets, sigma):
    """P1 stiffness matrix for div(sigma grad V) = 0, plus per-tet volumes."""
    p = nodes[tets] * 1e-3                       # mm -> m, so sigma*length is SI
    J = np.stack([p[:, 1] - p[:, 0], p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]], axis=2)
    det = np.linalg.det(J)
    vol = det / 6.0
    flip = vol < 0
    if flip.any():                               # enforce positive orientation
        tets = tets.copy()
        tets[flip] = tets[flip][:, [0, 2, 1, 3]]
        p = nodes[tets] * 1e-3
        J = np.stack([p[:, 1] - p[:, 0], p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]], axis=2)
        vol = np.linalg.det(J) / 6.0
    Jinv = np.linalg.inv(J)
    g = np.empty((len(tets), 4, 3))
    # x - p0 = J @ (lam1,lam2,lam3)  =>  lam = J^-1 (x - p0), so grad lam_i is
    # ROW i of J^-1, i.e. Jinv itself with no transpose. (Transposing here is a
    # silent error: it still gives a symmetric matrix with zero row sums, so
    # current balance still checks out while the field is wrong by ~2x. Caught
    # only by fem/scripts/verify_solver.py against the analytic sphere.)
    g[:, 1:, :] = Jinv
    g[:, 0, :] = -g[:, 1:, :].sum(axis=1)
    Ke = (sigma * vol)[:, None, None] * np.einsum("eik,ejk->eij", g, g)
    I = np.repeat(tets, 4, axis=1).ravel()
    Jj = np.tile(tets, (1, 4)).ravel()
    A = sp.coo_matrix((Ke.ravel(), (I, Jj)), shape=(len(nodes),) * 2).tocsr()
    return A, vol, tets


def boundary_nodes(tets):
    """Nodes on faces owned by exactly one tetrahedron, with their face areas."""
    f = np.concatenate([tets[:, [0, 1, 2]], tets[:, [0, 1, 3]],
                        tets[:, [0, 2, 3]], tets[:, [1, 2, 3]]])
    key = np.sort(f, axis=1)
    _, idx, cnt = np.unique(key, axis=0, return_index=True, return_counts=True)
    return f[idx[cnt == 1]]



def dura_leak_report(nodes, tets, lab):
    """How much of the dura barrier the mesh actually resolves.

    The dura is only ~0.5 mm thick and tissues are assigned per TETRAHEDRON, so
    wherever the local element is too coarse the dura layer can be punched
    through and epidural fat ends up face-to-face with CSF -- a short across the
    most resistive tissue in the model.  This measures that directly instead of
    hoping it did not happen: it reports the area of every internal face whose
    two tetrahedra are (epidural or lead) on one side and (CSF, white or grey)
    on the other, against the total area of the faces that do sit on a dura
    boundary.
    """
    idx = {n: i for i, n in enumerate(C.ORDER)}
    outer = {idx["epidural"], idx["insulator"]} | {idx["contact%d" % i] for i in range(1, 9)}
    inner = {idx["csf"], idx["white"], idx["grey"]}
    dura = idx["dura"]

    f = np.concatenate([tets[:, [0, 1, 2]], tets[:, [0, 1, 3]],
                        tets[:, [0, 2, 3]], tets[:, [1, 2, 3]]])
    owner = np.tile(np.arange(len(tets)), 4)
    key = np.sort(f, axis=1)
    o = np.lexsort((key[:, 2], key[:, 1], key[:, 0]))
    key, f, owner = key[o], f[o], owner[o]
    same = np.all(key[1:] == key[:-1], axis=1)
    i0 = np.flatnonzero(same)
    la, lb = lab[owner[i0]], lab[owner[i0 + 1]]
    p = nodes[f[i0]]
    area = 0.5 * np.linalg.norm(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1)

    ina = np.isin(la, list(inner))
    inb = np.isin(lb, list(inner))
    oua = np.isin(la, list(outer))
    oub = np.isin(lb, list(outer))
    leak = area[(ina & oub) | (inb & oua)].sum()
    dwall = area[((la == dura) & (inb | oub)) | ((lb == dura) & (ina | oua))].sum()
    print("  dura barrier: resolved interface %.1f mm2, LEAK (epidural|CSF direct) "
          "%.3f mm2 = %.3f %% of it" % (dwall, leak, 100 * leak / max(dwall, 1e-12)),
          flush=True)

def main():
    t0 = time.time()
    d = np.load(os.path.join(C.OUT, "mesh.npz"))
    nodes, tets = d["nodes"], d["tets"]
    print("mesh: %d nodes, %d tets" % (len(nodes), len(tets)), flush=True)

    lab = classify_tets(nodes, tets)
    print("classified %.0fs" % (time.time() - t0), flush=True)
    names = C.ORDER + ["background"]
    cnt = {names[i] if i >= 0 else "background": int((lab == i).sum())
           for i in list(range(len(C.ORDER))) + [-1]}
    print("tets by tissue:", cnt, flush=True)

    keep = np.ones(len(tets), dtype=bool) if INCLUDE_BACKGROUND else (lab >= 0)
    tets_k, lab_k = tets[keep], lab[keep]
    used = np.unique(tets_k)
    remap = np.full(len(nodes), -1, dtype=np.int64)
    remap[used] = np.arange(len(used))
    nodes_k = nodes[used]
    tets_k = remap[tets_k]
    print("active: %d nodes, %d tets (background %s)"
          % (len(nodes_k), len(tets_k), "IN" if INCLUDE_BACKGROUND else "REMOVED"), flush=True)

    dura_leak_report(nodes_k, tets_k, lab_k)

    sigma = sigma_of(lab_k)
    A, vol, tets_k = assemble(nodes_k, tets_k, sigma)
    print("assembled %.0fs  nnz=%d" % (time.time() - t0, A.nnz), flush=True)

    # nodal volume shares, used to spread the injected current inside a contact
    nodal = np.zeros(len(nodes_k))
    np.add.at(nodal, tets_k.ravel(), np.repeat(vol / 4.0, 4))

    # -1 A return, spread over the insulating outer boundary by nodal area
    bf = boundary_nodes(tets_k)
    pb = nodes_k[bf] * 1e-3
    area = 0.5 * np.linalg.norm(np.cross(pb[:, 1] - pb[:, 0], pb[:, 2] - pb[:, 0]), axis=1)
    nodal_area = np.zeros(len(nodes_k))
    np.add.at(nodal_area, bf.ravel(), np.repeat(area / 3.0, 3))
    ret = nodal_area / nodal_area.sum()

    srcs = {}
    for i in range(1, 9):
        m = lab_k == C.ORDER.index("contact%d" % i)
        w = np.zeros(len(nodes_k))
        np.add.at(w, tets_k[m].ravel(), np.repeat(vol[m] / 4.0, 4))
        if w.sum() <= 0:
            raise SystemExit("contact %d got no elements -- mesh too coarse" % i)
        srcs[i] = w / w.sum()

    # gauge: pin the node furthest from the lead axis
    axis = np.array([55.60, 87.0])
    far = int(np.argmax(np.linalg.norm(nodes_k[:, :2] - axis, axis=1)))
    free = np.ones(len(nodes_k), dtype=bool)
    free[far] = False
    fi = np.flatnonzero(free)
    Aff = A[fi][:, fi].tocsr()

    import pyamg
    ml = pyamg.smoothed_aggregation_solver(Aff, max_coarse=500)
    print("AMG built %.0fs\n%s" % (time.time() - t0, ml.levels[0].A.shape), flush=True)

    phi = np.zeros((8, len(nodes_k)))
    for i in range(1, 9):
        b = (srcs[i] - ret)[fi]
        res = []
        x = ml.solve(b, tol=1e-11, maxiter=400, accel="cg", residuals=res)
        v = np.zeros(len(nodes_k))
        v[fi] = x
        phi[i - 1] = v
        rel = res[-1] / max(res[0], 1e-300)
        react = float(np.asarray(A[far].dot(v)).ravel()[0] - (srcs[i][far] - ret[far]))
        print("  contact %d: %3d its, rel resid %.2e, pin reaction %.2e A"
              % (i, len(res) - 1, rel, react), flush=True)

    V = phi[C.SOURCE_CONTACT - 1] - phi[C.SINK_CONTACT - 1]
    V *= C.CURRENT_A

    # ---- verification, all measured not assumed ----
    print("\n--- checks ---")
    for i in range(1, 9):
        m = lab_k == C.ORDER.index("contact%d" % i)
        nd = np.unique(tets_k[m])
        spread = V[nd].max() - V[nd].min()
        role = {C.SOURCE_CONTACT: "SOURCE", C.SINK_CONTACT: "SINK"}.get(i, "float")
        print("  contact %d %-6s  V = %+9.4f V   spread %.3e V (%.4f %% of drive)"
              % (i, role, V[nd].mean(), spread, 100 * spread / (V.max() - V.min())))
    # net current through every contact: sum of A@V restricted to its interior nodes
    r = A.dot(V)
    for i in range(1, 9):
        m = lab_k == C.ORDER.index("contact%d" % i)
        nd = np.unique(tets_k[m])
        # interior-of-body nodes only touch that body, so sum(r) over body nodes
        # is the net current entering it
        print("  contact %d net current %+.6f A" % (i, r[nd].sum()))
    print("  global sum of nodal currents %.3e A (should be ~0)" % r.sum())
    E = np.abs(np.diff(np.sort(V)[[0, -1]]))
    print("  V range: %.4f V  (max %.4f, min %.4f)" % (E[0], V.max(), V.min()))

    np.savez_compressed(
        os.path.join(C.OUT, "solution_bg.npz" if INCLUDE_BACKGROUND else "solution.npz"),
        nodes=nodes_k, tets=tets_k, label=lab_k, sigma=sigma,
        phi=phi.astype(np.float32), V=V.astype(np.float32),
        order=np.array(C.ORDER), pinned_node=far,
        background_included=INCLUDE_BACKGROUND)
    print("wrote %s  %.0fs" % ("solution_bg.npz" if INCLUDE_BACKGROUND
                                 else "solution.npz", time.time() - t0))


if __name__ == "__main__":
    main()
