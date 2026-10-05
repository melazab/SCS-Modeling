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
from mesh_faces import face_batches
import artifacts

INCLUDE_BACKGROUND = os.environ.get("SCS_BACKGROUND", "0") == "1"


CHUNK = 400_000
ACCEPTED_RESIDUAL = 1e-6


def check_residual(matrix, x, rhs, contact, iterations):
    """Check the actual equation, independently of CG's recursive residual."""
    rel = float(np.linalg.norm(matrix @ x - rhs) / max(np.linalg.norm(rhs), 1e-300))
    if not np.isfinite(x).all() or not np.isfinite(rel) or rel > ACCEPTED_RESIDUAL:
        raise RuntimeError("contact %d failed convergence: true relative residual %.3e "
                           "exceeds accepted tolerance %.1e after %d iterations; "
                           "no solution published" % (contact, rel, ACCEPTED_RESIDUAL, iterations))
    return rel


def solve_contact(ml, matrix, rhs, contact):
    """Continue capped CG batches until accepted, with a firm 1600-step limit."""
    x = None
    iterations = 0
    for attempt in range(4):
        history = []
        x = ml.solve(rhs, x0=x, tol=1e-11, maxiter=400, accel='cg', residuals=history)
        iterations += max(0, len(history) - 1)
        try:
            rel = check_residual(matrix, x, rhs, contact, iterations)
            return x, rel, iterations
        except RuntimeError:
            if attempt == 3 or not np.isfinite(x).all():
                raise
            print('  contact %d: residual above acceptance after %d iterations; continuing'
                  % (contact, iterations), flush=True)


def build_testers(contact_stl=None, insulator_stl=None, order=None, progress=False):
    """{order name: [InsideTester, ...]} -- ONE TESTER PER STL BODY, grouped
    by the tissue name it belongs to.

    A tissue name used to mean exactly one shell (config.TISSUE_STL's five).
    Since 2026-09-15 it means a CLASS (config.TISSUE_BODIES): "root" is 123
    separate root/rootlet shells, "dura" is the main meninges plus 24 root
    sheaths and DRG coatings, and so on. A point is inside the class if it
    is inside ANY of its bodies.

    The bodies of a class are deliberately NOT concatenated into one
    InsideTester. InsideTester decides by +Z ray PARITY, so a point inside
    two overlapping shells of the same class would count two crossings and
    come back OUTSIDE. RADO's shells within a class do overlap (e.g. a root
    branch and the rootlet it splits into), so merging would silently punch
    holes in exactly the tissues this change is adding."""
    contact_stl = C.CONTACT_STL if contact_stl is None else contact_stl
    insulator_stl = C.INSULATOR_STL if insulator_stl is None else insulator_stl
    order = C.order_for(contact_stl) if order is None else order
    paths = {"contact%d" % i: [p] for i, p in contact_stl.items()}
    paths["insulator"] = [insulator_stl]
    for k, v in C.TISSUE_BODIES.items():
        paths[k] = list(v)
    testers = {}
    n_bodies = 0
    for name in order:
        testers[name] = [InsideTester(*read_stl(p)) for p in paths[name]]
        n_bodies += len(testers[name])
    if progress:
        print("  built %d point-in-shell testers over %d tissue names"
              % (n_bodies, len(order)), flush=True)
    return testers


def classify_tets(nodes, tets, contact_stl=None, insulator_stl=None, order=None,
                  testers=None):
    """Label every tet by the tissue containing its centroid (ORDER: first wins).

    `contact_stl`/`insulator_stl` override which lead's geometry the contact%d/
    insulator testers are built from -- default to `config.CONTACT_STL`/
    `config.INSULATOR_STL` (the frozen fem/leads/fem_dorsal_T10 lead) if
    omitted, so every existing call site (this module's own main(), and
    mesh_preview.py's default-lead preview) is unchanged. Tissue anatomy
    (config.TISSUE_BODIES) is never overridden -- only the lead is swappable.

    `order` is the classification priority list (see config.order_for());
    defaults to `config.order_for(contact_stl)`, i.e. exactly as many
    "contactN" entries as `contact_stl` actually has keys for -- so an
    arbitrary N-contact lead never gets matched against config.ORDER's
    fixed 8-contact block (which would KeyError on testers["contact5"] for
    a 4-contact lead, or silently misclassify tets for any other N).

    `testers` lets a caller that already built them (build_testers()) reuse
    them instead of re-reading 245 STLs; default None builds them here.

    PERFORMANCE, MEASURED 2026-09-15 on the 1 887 693-tet field-of-record
    mesh (fem/out/mesh.npz), same machine, back to back. This was the flagged
    risk of going from 14 testers to 249:

        5 tissue names,  14 bodies  (the old set)       6.7 s
       11 tissue names, 249 bodies, no bbox skip       27.6 s  (4.8 build
                                                               + 22.9 classify)
       11 tissue names, 249 bodies, with bbox skip     25.8 s  (4.6 build
                                                               + 21.2 classify)

    ~3.9x slower in absolute terms, and ~19 s added to a build whose size
    field alone runs 12-18 minutes -- so it is not the bottleneck and needs
    no mitigation. On the 66 939-tet coarse mesh it is 1.2 s -> 6.5 s, almost
    all of which is the one-off 4.6 s of reading 249 STLs and building their
    XY grids; pass `testers=` to reuse them across meshes.

    The bounding-box skip below is a small, exactly-equivalent win (22.9 ->
    21.2 s, identical labels verified by comparing the full label arrays):
    a body that cannot contain any point of this chunk is skipped before
    InsideTester's own 6-comparison inbox pass over the chunk. All 32
    sympathetic-chain bodies, for example, sit entirely outside the mesh box.
    """
    order = C.order_for(C.CONTACT_STL if contact_stl is None else contact_stl) \
        if order is None else order
    if testers is None:
        testers = build_testers(contact_stl, insulator_stl, order)
    cen = np.empty((len(tets), 3), dtype=np.float64)
    for start in range(0, len(tets), CHUNK):
        cen[start:start + CHUNK] = nodes[tets[start:start + CHUNK]].mean(axis=1)
    lab = np.full(len(cen), -1, dtype=np.int16)
    todo = np.arange(len(cen))
    for i, name in enumerate(order):
        if todo.size == 0:
            break
        hit = np.zeros(todo.size, dtype=bool)
        for s in range(0, todo.size, CHUNK):
            sl = slice(s, s + CHUNK)
            pts = cen[todo[sl]]
            plo, phi = pts.min(axis=0), pts.max(axis=0)
            sub = hit[sl]
            for t in testers[name]:
                if np.any(t.hi < plo - 1e-9) or np.any(t.lo > phi + 1e-9):
                    continue          # this body cannot contain any point here
                sub |= t(pts)
            hit[sl] = sub
        print("  classified %s: %d tets" % (name, hit.sum()), flush=True)
        lab[todo[hit]] = i
        todo = todo[~hit]
    return lab


def sigma_of(lab, include_background=None, order=None):
    """`order` is the SAME list classify_tets() used to produce `lab`
    (default config.ORDER, the frozen lead) -- must match, since `lab`'s
    integer values are positions into it. Already generic over N: a
    "contactN" name for ANY N gets SIGMA_METAL via the startswith check
    below, nothing here assumed exactly 8 of them."""
    include_background = INCLUDE_BACKGROUND if include_background is None else include_background
    order = C.ORDER if order is None else order
    s = np.full(len(lab), C.SIGMA["background"] if include_background else np.nan)
    for i, name in enumerate(order):
        v = C.SIGMA_METAL if name.startswith("contact") else C.SIGMA[name]
        s[lab == i] = v
    return s


def assemble(nodes, tets, sigma, batch=None):
    """Bounded-memory reference assembly; production solves use Elmer."""
    from element_batches import assemble as batched_assemble
    return batched_assemble(nodes, tets, sigma, batch=batch)


def boundary_nodes(tets):
    """Nodes on faces owned by exactly one tetrahedron, with their face areas."""
    faces = [f[lone] for f, owner, paired, lone in face_batches(tets)]
    return np.concatenate(faces) if faces else np.empty((0, 3), dtype=np.int64)



def dura_leak_report(nodes, tets, lab, order=None):
    """How much of the dura barrier the mesh actually resolves.

    The dura is only ~0.5 mm thick and tissues are assigned per TETRAHEDRON, so
    wherever the local element is too coarse the dura layer can be punched
    through and epidural fat ends up face-to-face with CSF -- a short across the
    most resistive tissue in the model.  This measures that directly instead of
    hoping it did not happen: it reports the area of every internal face whose
    two tetrahedra are (epidural or lead) on one side and (CSF, white or grey)
    on the other, against the total area of the faces that do sit on a dura
    boundary.

    `order` is the SAME list `lab` was classified against (default
    config.ORDER, the frozen lead) -- see classify_tets()'s own docstring.
    The "every contact is outer" set below is now discovered from `order`
    itself (any name starting with "contact"), not from a hardcoded
    range(1, 9), so this works unchanged for an N-contact lead.

    2026-09-15: the tissue set grew from five compartments to eleven, and
    this metric's definition is DELIBERATELY UNCHANGED -- `outer` is still
    exactly epidural + lead, `inner` still exactly CSF + white + grey. The
    new tissues (vertebra, disc, root, blood, drg, sympathetic) are in
    neither set. Two reasons: this is the number fem/README.md publishes and
    the Mesh Generator panel displays, so it has to stay comparable across
    that change; and the new tissues genuinely sit on BOTH sides of the dura
    (a root pierces it, a vessel runs through fat and through the
    subarachnoid space), so there is no honest way to put them on one side
    without inventing an answer. What this measures is still what it always
    measured: is the ~0.5 mm dura barrier resolved, or does epidural fat
    short straight to CSF. Measured after the change, on the field-of-record
    mesh: 0.000 mm2, 0.000 %, unchanged.
    """
    order = C.ORDER if order is None else order
    idx = {n: i for i, n in enumerate(order)}
    outer = {idx["epidural"], idx["insulator"]} | {v for k, v in idx.items() if k.startswith("contact")}
    inner = {idx["csf"], idx["white"], idx["grey"]}
    dura = idx["dura"]

    leak = dwall = 0.0
    for f, owner, i0, lone in face_batches(tets):
        la, lb = lab[owner[i0]], lab[owner[i0 + 1]]
        ina, inb = np.isin(la, list(inner)), np.isin(lb, list(inner))
        oua, oub = np.isin(la, list(outer)), np.isin(lb, list(outer))
        leaks = (ina & oub) | (inb & oua)
        walls = ((la == dura) & (inb | oub)) | ((lb == dura) & (ina | oua))
        selected = leaks | walls
        # Only material interfaces need coordinates; the old code expanded
        # EVERY interior face to a (n,3,3) float64 array and exhausted RAM.
        p = nodes[f[i0[selected]]]
        area = 0.5 * np.linalg.norm(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1)
        leak += float(area[leaks[selected]].sum())
        dwall += float(area[walls[selected]].sum())
    pct = 100 * leak / max(dwall, 1e-12)
    print("  dura barrier: resolved interface %.1f mm2, LEAK (epidural|CSF direct) "
          "%.3f mm2 = %.3f %% of it" % (dwall, leak, pct),
          flush=True)
    return dict(dura_wall_mm2=float(dwall), leak_mm2=float(leak), leak_pct=float(pct))

def solve_basis_fields(nodes, tets, contact_stl=None, insulator_stl=None,
                        include_background=None, t0=None, progress_cb=None):
    """Classify tissue, assemble div(sigma grad V) = 0, and solve N basis
    fields (N = the number of contacts in `contact_stl`, 8 for the frozen
    lead) for an arbitrary (nodes, tets) mesh and, optionally, an arbitrary
    lead's contact/insulator geometry -- e.g. a fresh mesh built by
    SCS_Mesh_Generator.FCMacro's Stage 1 against a Lead-Designer-exported
    lead of ANY contact count, rather than the frozen fem/leads/fem_dorsal_T10
    this module's own main() always used.

    This is the WHOLE physics body that used to live directly in main() --
    extracted so main() (fem/run_all.sh's field of record, and everyone who
    reads solution.npz) and the new interactive "Solve Field" path
    (fem/scripts/solve_lead.py) share the exact same code, not two copies
    that can drift. main() below calls this with every argument at its old
    default (`contact_stl`/`insulator_stl` = None = config.py's frozen
    8-contact lead, `include_background` = None = module INCLUDE_BACKGROUND),
    so its result and every number it prints are UNCHANGED from before this
    refactor AND from before this function's own generalization to an
    arbitrary contact count -- verified by re-running against
    fem/out/mesh.npz and diffing the resulting phi/V against the committed
    fem/out/solution.npz (see fem/scripts/solve_lead.py's module docstring
    and this session's report).

    Returns a dict with exactly the fields main() has always written to
    solution.npz (nodes, tets, label, sigma, phi, V, order, pinned_node,
    background_included) -- np.savez_compressed(path, **result) reproduces
    the old call exactly. `phi` is shaped (N, nnode) and `order` is now
    N-long's own contact block (config.order_for(contact_stl)), not always
    8 -- both were hardcoded to 8 before this generalization.

    `V` (the SOURCE_CONTACT/SINK_CONTACT bipolar reference field) is only
    computed when this lead actually HAS both of those contact numbers
    (true for the frozen 8-contact lead: SOURCE=3, SINK=5). For a lead with
    fewer contacts than SINK_CONTACT, `V` is a documented zero placeholder
    instead -- config.SOURCE_CONTACT/SINK_CONTACT are meaningful only for
    the frozen lead's own field-of-record case, not for an arbitrary live
    lead, which gets its currents from scs_montage_feature.py's per-contact
    properties applied directly to `phi`, never from V.

    `progress_cb`, if given, is called with an int 0-100 at a few
    milestones: classification done, assembly done, AMG setup done, then
    once per basis-contact solve (the N CG solves are roughly even in cost,
    so these are evenly spaced) -- per Mohamed's brief ("AMG built, then
    after each of the 8 basis solves for the solve stage" -- "8" there was
    the frozen lead being discussed at the time; the spacing itself was
    always N-based, see the `_progress` call inside the solve loop below).
    Default None (main() below passes nothing) so this is purely additive;
    solve_lead.py is the one caller that passes a real callback.
    """
    def _progress(pct):
        if progress_cb is not None:
            progress_cb(pct)

    t0 = time.time() if t0 is None else t0
    include_background = INCLUDE_BACKGROUND if include_background is None else include_background
    contact_stl = C.CONTACT_STL if contact_stl is None else contact_stl
    insulator_stl = C.INSULATOR_STL if insulator_stl is None else insulator_stl
    order = C.order_for(contact_stl)
    contact_ids = sorted(contact_stl)

    import mesh_cost
    available = mesh_cost.available_ram_gb()
    estimated_gb = 1.0 + len(tets) * 2500 / 1e9
    if available is not None and estimated_gb > 0.8 * available:
        raise RuntimeError('Solve needs approximately %.1f GB; %.1f GB available. '
                           'Use a coarser mesh or a larger machine.' % (estimated_gb, available))

    lab = classify_tets(nodes, tets, contact_stl=contact_stl,
                         insulator_stl=insulator_stl, order=order)
    print("classified %.0fs" % (time.time() - t0), flush=True)
    _progress(5)
    names = order + ["background"]
    cnt = {names[i] if i >= 0 else "background": int((lab == i).sum())
           for i in list(range(len(order))) + [-1]}
    print("tets by tissue:", cnt, flush=True)

    keep = np.ones(len(tets), dtype=bool) if include_background else (lab >= 0)
    tets_k, lab_k = tets[keep], lab[keep]
    used = np.unique(tets_k)
    remap = np.full(len(nodes), -1, dtype=np.int64)
    remap[used] = np.arange(len(used))
    nodes_k = nodes[used]
    tets_k = remap[tets_k]
    print("active: %d nodes, %d tets (background %s)"
          % (len(nodes_k), len(tets_k), "IN" if include_background else "REMOVED"), flush=True)

    dura_leak_report(nodes_k, tets_k, lab_k, order=order)

    sigma = sigma_of(lab_k, include_background=include_background, order=order)
    A, vol, tets_k = assemble(nodes_k, tets_k, sigma)
    from scipy.sparse.csgraph import connected_components
    n_components = connected_components(A, directed=False, return_labels=False)
    if n_components != 1:
        raise RuntimeError("active mesh has %d disconnected components; a single gauge "
                           "cannot define this solve. Refine/check tissue connectivity." % n_components)
    print("assembled %.0fs  nnz=%d" % (time.time() - t0, A.nnz), flush=True)
    _progress(8)

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
    for i in contact_ids:
        m = lab_k == order.index("contact%d" % i)
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
    np.random.seed(0)  # reproducible AMG spectral estimates / setup
    ml = pyamg.smoothed_aggregation_solver(Aff, max_coarse=500)
    print("AMG built %.0fs\n%s" % (time.time() - t0, ml.levels[0].A.shape), flush=True)
    _progress(12)

    n_contacts = len(contact_ids)
    residuals = []
    iteration_counts = []
    phi = np.zeros((n_contacts, len(nodes_k)))
    for pos, i in enumerate(contact_ids):
        b = (srcs[i] - ret)[fi]
        x, rel, iterations = solve_contact(ml, Aff, b, i)
        v = np.zeros(len(nodes_k))
        v[fi] = x
        phi[pos] = v
        residuals.append(rel)
        iteration_counts.append(iterations)
        react = float(np.asarray(A[far].dot(v)).ravel()[0] - (srcs[i][far] - ret[far]))
        print("  contact %d: %3d its, rel resid %.2e, pin reaction %.2e A"
              % (i, iterations, rel, react), flush=True)
        _progress(round(12 + 88 * (pos + 1) / n_contacts))

    if C.SOURCE_CONTACT in contact_ids and C.SINK_CONTACT in contact_ids:
        V = phi[contact_ids.index(C.SOURCE_CONTACT)] - phi[contact_ids.index(C.SINK_CONTACT)]
        V *= C.CURRENT_A
    else:
        # This lead doesn't have both of config.SOURCE_CONTACT/SINK_CONTACT
        # (meaningful only for the frozen 8-contact lead) -- V is a
        # documented zero placeholder, not a bipolar reference field. Every
        # real per-contact-current combination for THIS lead goes through
        # `phi` directly (scs_montage.py), never through V.
        V = np.zeros(len(nodes_k), dtype=np.float64)
        print("  NOTE: this lead has no contact %d and/or %d "
              "(config.SOURCE_CONTACT/SINK_CONTACT) among its %d contact(s) "
              "%s -- 'V' is a zero placeholder, not a bipolar reference "
              "field; use 'phi' with scs_montage.py for this lead's actual "
              "per-contact-current fields."
              % (C.SOURCE_CONTACT, C.SINK_CONTACT, n_contacts, contact_ids), flush=True)

    # ---- verification, all measured not assumed ----
    print("\n--- checks ---")
    v_drive = max(V.max() - V.min(), 1e-300)
    for i in contact_ids:
        m = lab_k == order.index("contact%d" % i)
        nd = np.unique(tets_k[m])
        spread = V[nd].max() - V[nd].min()
        role = {C.SOURCE_CONTACT: "SOURCE", C.SINK_CONTACT: "SINK"}.get(i, "float")
        print("  contact %d %-6s  V = %+9.4f V   spread %.3e V (%.4f %% of drive)"
              % (i, role, V[nd].mean(), spread, 100 * spread / v_drive))
    # net current through every contact: sum of A@V restricted to its interior nodes
    r = A.dot(V)
    for i in contact_ids:
        m = lab_k == order.index("contact%d" % i)
        nd = np.unique(tets_k[m])
        # interior-of-body nodes only touch that body, so sum(r) over body nodes
        # is the net current entering it
        print("  contact %d net current %+.6f A" % (i, r[nd].sum()))
    print("  global sum of nodal currents %.3e A (should be ~0)" % r.sum())
    E = np.abs(np.diff(np.sort(V)[[0, -1]]))
    print("  V range: %.4f V  (max %.4f, min %.4f)" % (E[0], V.max(), V.min()))

    return dict(
        nodes=nodes_k, tets=tets_k, label=lab_k, sigma=sigma,
        phi=phi, V=V,
        relative_residuals=np.asarray(residuals), accepted_residual=ACCEPTED_RESIDUAL,
        iterations=np.asarray(iteration_counts), amg_seed=0,
        requested_residual=1e-11, contact_ids=np.asarray(contact_ids),
        order=np.array(order), pinned_node=far,
        background_included=include_background)


def main():
    t0 = time.time()
    d = np.load(os.path.join(C.OUT, "mesh.npz"))
    nodes, tets = d["nodes"], d["tets"]
    print("mesh: %d nodes, %d tets" % (len(nodes), len(tets)), flush=True)

    result = solve_basis_fields(nodes, tets, include_background=INCLUDE_BACKGROUND, t0=t0)

    out_name = "solution_bg.npz" if INCLUDE_BACKGROUND else "solution.npz"
    artifacts.atomic_npz(os.path.join(C.OUT, out_name), **result)
    print("wrote %s  %.0fs" % (out_name, time.time() - t0))


if __name__ == "__main__":
    main()
