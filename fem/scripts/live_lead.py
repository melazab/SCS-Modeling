"""Find and tessellate the live lead, preserving world-space geometry.

Run keys combine lead geometry, mesh parameters and model/code provenance.
Each run uses its own lead STL directory so another dock cannot replace inputs
mid-build. Only completed, provenance-checked artifacts are cache hits.
"""
import glob
import hashlib
import json
import os
import re
import struct

import numpy as np

import config as C

# 2026-09-14 CORRECTION: this used to be a hardcoded constant (N_CONTACTS =
# 8) that find_live_lead() enforced as the ONLY acceptable contact count,
# and tessellate_lead_arrays()/fingerprint_from_arrays()/write_stls() all
# looped range(1, N_CONTACTS + 1) instead of over whatever contacts the
# live lead actually had. That is exactly the bug Mohamed hit: he built a
# real 4-contact lead in Lead Designer (whose own "Contacts" spinbox
# already supports an arbitrary count) and got "Lead 1 has 4/8 contacts --
# the FEM pipeline only supports an 8-contact lead" -- a limitation of
# THIS file, not of Lead Designer or of the FEM physics (config.ORDER/
# assign_and_solve.solve_basis_fields()/scs_montage.py are now all generic
# over N too; see their own docstrings). Every function below now derives N
# from the live lead's OWN discovered contacts (find_live_lead()'s
# `n_contacts`/`contacts` dict), never from a fixed module constant. The
# name is kept as a documented reference to the FROZEN field-of-record
# lead's own contact count (fem/leads/fem_dorsal_T10 -- e.g. a caller that
# wants "how many contacts does the frozen lead have" for a message or a
# default), NOT as a requirement any live lead must match.
FROZEN_LEAD_N_CONTACTS = 8

# Comfortably inside the measured 0.02-0.3mm plateau that reproduces
# fem/leads/fem_dorsal_T10's own STL facet counts exactly -- see module
# docstring "TESSELLATION TOLERANCE" for the actual measurements.
TESSELLATE_TOL_MM = 0.1

LIVE_LEAD_DIR = os.path.join(C.LEAD_RUNS, "_live")

_CONTACT_RE = re.compile(r"^SCS_Preview_Lead(\d+)_Contact_(\d+)$")
_INSULATOR_RE = re.compile(r"^SCS_Preview_Lead(\d+)_Insulator$")


class LiveLeadStatus:
    OK = "ok"
    NONE = "none"
    MULTIPLE = "multiple"
    WRONG_SHAPE = "wrong_shape"


def find_live_lead(doc):
    """Discover every Lead Designer lead, including hidden objects.

    Solver IDs are consecutive in sorted (lead number, contact number) order.
    contact_map preserves the physical identity; all insulation shells share
    one material surface without joining their disconnected geometries.
    """
    if doc is None:
        return dict(status=LiveLeadStatus.NONE,
                    reason="No document open -- place a lead with Lead Designer first.")

    contacts_by_lead = {}
    insulators_by_lead = {}
    for obj in doc.Objects:
        m = _CONTACT_RE.match(obj.Name)
        if m:
            n, idx = int(m.group(1)), int(m.group(2))
            contacts_by_lead.setdefault(n, {})[idx] = obj
            continue
        m = _INSULATOR_RE.match(obj.Name)
        if m:
            insulators_by_lead[int(m.group(1))] = obj

    lead_indices = sorted(set(contacts_by_lead) | set(insulators_by_lead))
    if not lead_indices:
        return dict(status=LiveLeadStatus.NONE,
                    reason="No lead found -- place one with Lead Designer first.")
    contacts, contact_map = {}, {}
    for n in lead_indices:
        have = sorted(contacts_by_lead.get(n, {}))
        if not have or have != list(range(1, len(have)+1)) or n not in insulators_by_lead:
            return dict(status=LiveLeadStatus.WRONG_SHAPE,
                        reason="Lead %d needs consecutive contacts starting at 1 and an insulator." % n)
        for local in have:
            i = len(contacts)+1
            contacts[i] = contacts_by_lead[n][local]
            contact_map[i] = [n, local]
    return dict(status=LiveLeadStatus.OK, reason="", lead_index=lead_indices[0],
                lead_indices=lead_indices, contacts=contacts, contact_map=contact_map,
                insulators=insulators_by_lead,
                insulator=insulators_by_lead[lead_indices[0]], n_contacts=len(contacts))


def ready_text(info):
    return "%d lead(s) ready (%d contacts + %d insulators)." % (
        len(info['lead_indices']), info['n_contacts'], len(info['insulators']))


def contact_label(info, i):
    lead, contact = info['contact_map'][i]
    return "Lead %d / contact %d" % (lead, contact)


class LeadArrays(dict):
    """Geometry plus the persisted physical-to-solver contact mapping."""
    pass


def _tessellate_obj(obj, tol):
    """(verts (N,3) float64, tris (M,3) int64), DOCUMENT GLOBAL coordinates
    -- see module docstring "PLACEMENT IS ALREADY BAKED IN" for why no
    separate Placement multiply happens here.

    TESSELLATES A COPY OF THE SHAPE, NEVER obj.Shape DIRECTLY -- found live
    2026-09-14, not assumed: FreeCAD caches a mesh on a TopoShape's own BRep
    representation, keyed to the FINEST deflection ever requested against
    THAT shape object (including FreeCAD's own viewport auto-tessellation
    the very first time the object was made visible, at whatever "Deviation"
    the 3D view uses -- nothing this module controls). A LATER call with a
    coarser-or-equal `tol` silently returns the CACHED mesh instead of
    honouring the tolerance just passed in -- verified against this exact
    live document: Contact_01 (which an earlier debugging step had forced
    to a fine 0.001mm deflection) still returned 5024 triangles when asked
    for tol=0.1, while Contacts 2-8 (never touched, only ever
    viewport-tessellated) returned a stale 100 triangles each at the SAME
    tol=0.1 -- neither number is what tol=0.1 actually produces (measured:
    500). `obj.Shape.copy()` is a fresh TopoShape with no such cache, so
    tessellating the copy is the only way this function's `tol` argument
    is ever actually honoured."""
    shape = getattr(obj, "Shape", None)
    if shape is not None and not shape.isNull():
        shape = shape.copy()
    if shape is None or shape.isNull():
        # A DRG lead previews as Mesh::Feature (build_lead_config.
        # preview_set(), "A copied DRG lead previews as Mesh::Feature");
        # not reachable for a dorsal/ventral lead's 8-contact layout that
        # find_live_lead() requires, but handled rather than crashing if
        # someone ever mixes the two.
        mesh = getattr(obj, "Mesh", None)
        if mesh is None:
            raise ValueError("%s has neither a usable .Shape nor a .Mesh" % obj.Name)
        pts = np.array([[p.x, p.y, p.z] for p in mesh.Points], dtype=np.float64)
        tris = np.array([list(f.PointIndices) for f in mesh.Facets], dtype=np.int64)
        return pts, tris
    verts, tris = shape.tessellate(tol)
    v = np.array([[p.x, p.y, p.z] for p in verts], dtype=np.float64)
    t = np.array(tris, dtype=np.int64)
    return v, t


def tessellate_lead_arrays(lead_info, tol=TESSELLATE_TOL_MM):
    """{1..N: (verts, tris), "insulator": (verts, tris)} for an OK
    lead_info (see find_live_lead()) -- N = lead_info["n_contacts"],
    whatever contact count THIS live lead actually has, not a fixed 8.
    Pure in-memory tessellation -- no file I/O, so this is cheap enough to
    call just to compute a fingerprint without committing to a mesh/solve
    run."""
    if lead_info["status"] != LiveLeadStatus.OK:
        raise ValueError("tessellate_lead_arrays() called on a non-OK "
                          "lead_info: %s" % lead_info["reason"])
    out = LeadArrays()
    out.contact_map = lead_info["contact_map"]
    for i in range(1, lead_info["n_contacts"] + 1):
        out[i] = _tessellate_obj(lead_info["contacts"][i], tol)
    vertices, triangles, offset = [], [], 0
    out.insulators = {}
    for lead, obj in lead_info['insulators'].items():
        v, t = _tessellate_obj(obj, tol)
        out.insulators[lead] = (v, t)
        vertices.append(v)
        triangles.append(t + offset)
        offset += len(v)
    out['insulator'] = (np.concatenate(vertices), np.concatenate(triangles))
    return out


def fingerprint_from_arrays(arrays):
    """sha256 hex digest (16 hex chars) over the ACTUAL tessellated
    geometry -- see module docstring "FINGERPRINTING" for why this, not
    Shape.hashCode(). Iterates over however many contact keys `arrays`
    actually has (sorted, so key order is deterministic), not a fixed
    range(1, 9) -- a DIFFERENT contact count always changes the set of
    keys hashed, so switching leads (4-contact to 6-contact, say) always
    changes the fingerprint even before any coordinate is compared."""
    h = hashlib.sha256()
    mapping = getattr(arrays, 'contact_map', {})
    if len({v[0] for v in mapping.values()}) > 1:
        h.update(json.dumps(mapping, sort_keys=True).encode())
    contact_ids = sorted(k for k in arrays if k != "insulator")
    for key in contact_ids + ["insulator"]:
        v, t = arrays[key]
        h.update(str(key).encode())
        h.update(np.round(v, 6).astype(np.float64).tobytes())
        h.update(np.ascontiguousarray(t, dtype=np.int64).tobytes())
    return h.hexdigest()[:16]


def run_key_for(geom_fingerprint, mesh_params):
    """The cache key for fem/out/lead_runs/<run_key>/: a function of BOTH
    the lead's geometry AND the mesh size-field parameters, so Mesh
    Generator experimenting with a finer H_Min on the SAME lead never
    reuses (or clobbers) a coarser mesh Potential Visualizer already built
    for that lead with build_mesh.DEFAULT_PARAMS, while an UNCHANGED
    params dict on the same geometry -- the common case -- lands in the
    exact same directory either tool built it from, which is the whole
    point of sharing one cache."""
    import artifacts
    return artifacts.digest(dict(geometry=geom_fingerprint, params=mesh_params,
                                 model=artifacts.model_signature()))[:24]


def _write_binary_stl(path, verts, tris):
    """Minimal binary STL writer -- verts (N,3), tris (M,3) int. Normals
    are recomputed from the winding (stlio.read_stl(), the only reader
    this pipeline uses, ignores stored normals entirely and only reads
    vertex positions, so an all-zero-normal facet is harmless if the
    winding is degenerate)."""
    tri_pts = verts[np.ascontiguousarray(tris, dtype=np.int64)].astype(np.float32)
    normals = np.cross(tri_pts[:, 1] - tri_pts[:, 0], tri_pts[:, 2] - tri_pts[:, 0])
    ln = np.linalg.norm(normals, axis=1, keepdims=True)
    normals = np.divide(normals, ln, out=np.zeros_like(normals), where=ln > 0).astype(np.float32)
    rec = np.zeros((len(tri_pts), 50), dtype=np.uint8)
    rec[:, 0:12] = normals.view(np.uint8).reshape(-1, 12)
    rec[:, 12:48] = tri_pts.reshape(len(tri_pts), 9).view(np.uint8).reshape(-1, 36)
    with open(path, "wb") as fh:
        fh.write(b"SCS live lead -- fem/scripts/live_lead.py, do not edit".ljust(80, b"\0")[:80])
        fh.write(struct.pack("<I", len(tri_pts)))
        fh.write(rec.tobytes())


def write_stls(arrays, out_dir=LIVE_LEAD_DIR):
    """Write `arrays` (tessellate_lead_arrays()' return -- N contacts,
    whatever N the live lead had) to `out_dir` using config.lead_stls()'
    own naming convention -- so out_dir is immediately usable as a
    `lead_dir` for config.lead_stls()/is_single_lead_dir(),
    build_mesh.build(contact_stl=, insulator_stl=), assign_and_solve.
    classify_tets()/solve_basis_fields(), and solve_lead.py's --lead-dir,
    with NO changes to any of that code -- the whole point of keeping the
    same directory-of-(N+1)-STLs contract. Returns out_dir.

    Clears any PRE-EXISTING "SCS Lead Electrode *.stl" files in `out_dir`
    first: this directory is documented (module docstring, CACHE LAYOUT) as
    fully overwritten on every click, but that was only actually true when
    every lead had the same N=8 -- switching from an 8-contact lead to a
    4-contact one used to leave contacts 5-8's stale STLs sitting there,
    where config.lead_stls()'s auto-detection (detect_contact_count()) or a
    stale directory listing would still pick them up as if they belonged to
    the new, smaller lead. Only files matching this lead's own naming
    pattern are touched -- never a directory this function didn't itself
    populate."""
    os.makedirs(out_dir, exist_ok=True)
    for stale in glob.glob(os.path.join(out_dir, "SCS Lead Electrode *.stl")):
        os.remove(stale)
    contact_ids = sorted(k for k in arrays if k != "insulator")
    contact_stl, insulator_stl = C.lead_stls(out_dir, n_contacts=len(contact_ids))
    for i in contact_ids:
        v, t = arrays[i]
        _write_binary_stl(contact_stl[i], v, t)
    v, t = arrays["insulator"]
    _write_binary_stl(insulator_stl, v, t)
    write_insulator_surfaces(arrays, out_dir, overwrite=True)
    import artifacts
    artifacts.atomic_json(os.path.join(out_dir, 'contact_map.json'),
                          getattr(arrays, 'contact_map', {}))
    return out_dir


def write_insulator_surfaces(arrays, lead_dir, overwrite=False):
    """Display-only per-lead surfaces; the combined solver insulation is unchanged."""
    for lead, (vertices, triangles) in getattr(arrays, 'insulators', {}).items():
        path = os.path.join(lead_dir, 'insulator_lead_%d.stl' % lead)
        if overwrite or not os.path.exists(path):
            _write_binary_stl(path, vertices, triangles)


def save_report(run_dir, name, report):
    """Persist a mesh_preview.py/solve_lead.py JSON report dict as
    <run_dir>/<name>.json, so a later cache hit (same run_key, no
    subprocess needed) can repaint the UI from the SAME numbers a real run
    produced, not placeholders. `name` is "mesh_report" or "solve_report"."""
    import artifacts
    artifacts.atomic_json(os.path.join(run_dir, name + '.json'), report)


def load_report(run_dir, name):
    """The report saved by save_report(), or None if this run_dir predates
    this caching mechanism (e.g. a mesh.npz written before this session) --
    callers fall back to a generic 'cached' status rather than crashing on
    a missing file that isn't actually required for correctness (the .npz
    files are what downstream code reads)."""
    path = os.path.join(run_dir, name + ".json")
    if not os.path.isfile(path):
        return None
    with open(path) as fh:
        return json.load(fh)
