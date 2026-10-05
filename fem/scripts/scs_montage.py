"""Arbitrary per-contact-current montage -> potential field, sampled onto the
real anatomical tissue surfaces (epidural space, dura, CSF, white matter, grey matter), for the
"SCS Potential Field Visualizer".

WHY THIS EXISTS

`assign_and_solve.py` already solved 8 basis fields into `fem/out/solution.npz`
-- `phi[i-1]` is the potential (V/A, transfer impedance) for +1 A into contact
i and -1 A spread over the insulating outer boundary (see fem/README.md
"Basis fields, for the NEURON stage"). Any zero-net-current montage is
therefore a pure linear combination,

    V = sum_i I_i(A) * phi[i-1]

with no new mesh and no new solve -- see field.py's own docstring, which
already states this contract. This module does exactly that arithmetic and
samples the result at the actual STL vertices of the selected "clinically
legible" tissue shells named in config.TISSUE_STL (epidural, dura, csf, white, grey),
rather than at the synthetic gmsh box the FEM ran on. It has NO FreeCAD
dependency -- it is importable and testable from a plain terminal (see
`__main__` below) -- and hands the FreeCAD-glue half
(scs_montage_feature.py / scs_potential_viz.FCMacro) two things per tissue:
a triangulated surface (verts+tris, identical to the source STL's own
triangulation) and one sampled V per vertex, written to a small .vtp file.

CORRECTNESS CHECK

`fem/out/solution.npz["V"]` is *already* phi[2]-phi[4] at 1 A each -- the
bipolar case documented in fem/README.md RESULTS (contact 3 +1297.67 V,
contact 5 -1369.91 V; that table's potentials are literally this array's
own extrema). `self_check()` below asks this module's *general* montage
path for that exact case and diffs it against the stored array node-for-
node, not just eyeballing a min/max. This is the one physics fact this
module must not get wrong -- see HANDOFF.md's own warning that a
conservation check ("current in equals current out") once passed a field
that was wrong by 2.2x; this checks the actual field, not a proxy for it.

COLOUR-LEGEND CLAMPING -- why every .vtp carries TWO point-data arrays

FreeCAD's Fem::FemPostPipeline / FemPostClipFilter ViewObject has no
Python-reachable min/max clamp for its colour legend. Checked directly
against FreeCAD's own C++ source (this build: ~/Repos/freecad,
Mod/Fem/Gui/ViewProviderFemPostObject.h/.cpp): the only paint-relevant
properties are Field, Component, Transparency, PlainColorEdgeOnSurface,
EdgeColor, NoneFieldColor, LineWidth, PointSize -- there is no Range/Min/Max
property. The legend's actual min/max is set by the *protected* C++ method
`ViewProviderFemPostObject::setRangeOfColorBar()`, called only from
`WriteColorData()` off the selected array's own data range (or the
-0.5/+0.5 default before any real data has been through it once -- the
"stale colour bar" bug this repo's other worker already hit and fixed by
toggling Field via "None"). `setRangeOfColorBar` is never exposed to
Python, and the underlying `Gui::SoFCColorBar::setRange()` lives on a
Coin3D scene node that isn't reachable from a document object at all. So
the legend cannot be told "always show -100..100" directly.

The workaround here is not a display trick but a second *data* array:
alongside the true "V_volts", every .vtp also carries "V_clamped_<N>V", the
same values hard-saturated to +/-N volts before they ever reach FreeCAD.
FreeCAD's autoscale then legitimately spans exactly +/-N, because that IS
the array's true range -- giving bulk-tissue gradient without lying about
what "V_volts" itself contains. Switching ViewObject.Field back to
"V_volts" recovers the true, unclamped number at any point; nothing is
thrown away, only what the legend autoscales *against* changes. N defaults
to the 95th percentile of |V| pooled over the three target tissues'
*own* vertices for the requested montage (see `choose_clamp_v`) -- checked
against this repo's actual data rather than assumed, per the brief.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

try:
    import scipy  # noqa: F401  (field.py needs scipy.spatial.cKDTree)
except ImportError:
    # FreeCAD's embedded Python (this repo's build: 3.14.4, matching
    # fem/.venv's own interpreter exactly) already has numpy and vtk on its
    # path via /usr/lib/python3/dist-packages, but not scipy/pyamg/gmsh --
    # those only exist in fem/.venv (see fem/README.md "Setup": the system
    # python is otherwise too new for some of these wheels). Since the venv
    # was built with --system-site-packages and the interpreter ABI matches
    # exactly, just adding its site-packages onto FreeCAD's sys.path is
    # sufficient -- no subprocess, no separate interpreter needed. If the
    # ABI ever stops matching, this import will fail loudly here instead of
    # silently picking up incompatible binaries.
    _venv_site = os.path.join(
        _HERE, "..", ".venv", "lib",
        "python%d.%d" % sys.version_info[:2], "site-packages",
    )
    sys.path.insert(0, os.path.abspath(_venv_site))
    import scipy  # noqa: F401  -- let this raise if it still fails

import numpy as np
import config as C
import stlio
from field import TetField

# 2026-09-14 CORRECTION: this used to be a hardcoded N_CONTACTS = 8 module
# constant, used to validate every currents_mA array's shape and to size
# default_bipolar_mA()'s array -- both wrong for a live lead with a
# different contact count. The actual contact count for ANY solution.npz is
# simply how many basis fields it has (phi.shape[0] / TetField.values.
# shape[0]), discovered per-call below rather than assumed fixed. Kept as a
# documented constant for the ONE case that genuinely is fixed: the frozen
# field-of-record lead (fem/leads/fem_dorsal_T10) always has 8 contacts,
# which is what default_bipolar_mA()'s own SOURCE/SINK convention (config.
# SOURCE_CONTACT=3, SINK_CONTACT=5) is defined against.
FROZEN_LEAD_N_CONTACTS = 8
CORE_TISSUES = ("epidural", "dura", "csf", "white", "grey")
TISSUE_LABELS = dict(epidural="Epidural space", dura="Dura", csf="CSF",
    white="White matter", grey="Grey matter", vertebra="Vertebrae", disc="Discs",
    blood="Vasculature", root="Nerve roots", drg="DRG", sympathetic="Sympathetic chain")
TISSUES = CORE_TISSUES + tuple(t for t in C.TISSUE_BODIES if t not in CORE_TISSUES)


def read_surface(paths):
    """Combine every anatomical body of a tissue class, preserving triangle indices."""
    if isinstance(paths, (str, os.PathLike)):
        return stlio.read_stl(paths)
    vertices, triangles, offset = [], [], 0
    for path in paths:
        v, tri = stlio.read_stl(path)
        vertices.append(v)
        triangles.append(tri + offset)
        offset += len(v)
    return np.concatenate(vertices), np.concatenate(triangles)


# Boundary points closer than this to the nearest solved FEM node get their
# value filled from that node instead of being left NaN -- see the long
# comment on FILL_TOL below and in sample_tissue(). 0.5 mm is the coarsest
# target element size anywhere in the cord/dura mesh (fem/README.md "The
# mesh"), so this never reaches past one element's width.
FILL_TOL_MM = 0.5

_CACHE = {}  # {(path, mtime): (TetField, npz_dict, node_cKDTree)} -- one entry


def release_unused_basis(active_paths):
    """Release only in-memory solutions not referenced by an open document."""
    active = {os.path.realpath(p) for p in active_paths if p}
    removed = 0
    for key in list(_CACHE):
        if os.path.realpath(key[0]) not in active:
            _field, archive, _tree = _CACHE.pop(key)
            archive.close()
            removed += 1
    return removed


def _load_basis(sol_npz=None):
    """Load solution.npz once per process and cache the TetField (its cKDTree
    over 1.7M tets is the expensive part -- ~seconds -- everything downstream,
    including repeated interactive montage edits, reuses this) plus a second,
    cheap cKDTree over the mesh's own 289136 nodes, used only to fill in
    boundary vertices TetField can't locate (see sample_tissue)."""
    sol_npz = sol_npz or os.path.join(C.OUT, "solution.npz")
    key = (sol_npz, os.path.getmtime(sol_npz))
    if key not in _CACHE:
        from scipy.spatial import cKDTree

        d = np.load(sol_npz, allow_pickle=True)
        f = TetField(d["nodes"], d["tets"], d["phi"].astype(np.float64))
        node_tree = cKDTree(d["nodes"])
        _CACHE.clear()
        _CACHE[key] = (f, d, node_tree)
    return _CACHE[key]


def default_bipolar_mA(n_contacts=FROZEN_LEAD_N_CONTACTS):
    """The existing validated case from config.py: SOURCE_CONTACT +CURRENT_A,
    SINK_CONTACT -CURRENT_A, everything else open. Same sign convention as
    assign_and_solve.py: positive = source/anode, negative = sink/cathode.

    `n_contacts` defaults to the frozen lead's own 8 -- config.SOURCE_CONTACT/
    SINK_CONTACT (3/5) are only meaningful for a lead that actually has
    those two contacts (see assign_and_solve.solve_basis_fields()'s own
    docstring); a caller with a different lead's solution.npz in hand
    should pass that lead's real contact count (e.g. its own phi.shape[0])
    so this raises a clear IndexError instead of silently building an
    8-long array for a 4-contact lead."""
    currents_mA = np.zeros(n_contacts)
    currents_mA[C.SOURCE_CONTACT - 1] = C.CURRENT_A * 1000.0
    currents_mA[C.SINK_CONTACT - 1] = -C.CURRENT_A * 1000.0
    return currents_mA


def currents_mA_to_A(currents_mA, n_contacts):
    """`n_contacts` is the ACTUAL number of basis fields the target
    solution.npz has (its phi.shape[0]) -- passed in by the caller (see
    montage_at_nodes/montage_at_points below), not assumed to be 8, so this
    validates correctly for a live N-contact lead's own solved field."""
    currents_mA = np.asarray(currents_mA, dtype=np.float64)
    if currents_mA.shape != (n_contacts,):
        raise ValueError(
            "expected %d per-contact currents (mA) for this lead, got shape %s"
            % (n_contacts, currents_mA.shape)
        )
    imbalance = float(currents_mA.sum())
    scale = max(1.0, float(np.abs(currents_mA).max()))
    if abs(imbalance) > 1e-6 * scale:
        print(
            "WARNING: montage currents do not sum to ~0 (sum = %.6g mA). "
            "This model is DC with no return path except through tissue -- "
            "an unbalanced montage is unphysical. Proceeding anyway "
            "(warning only, not blocking)." % imbalance
        )
    return currents_mA / 1000.0


def montage_at_nodes(currents_mA, sol_npz=None):
    """V (volts) at every one of the FEM mesh's own nodes, for arbitrary
    per-contact currents (as many as this solution.npz has basis fields
    for -- 8 for the frozen lead, N for a live N-contact one). Used by
    self_check() (frozen-lead consistency check) and by the live-lead
    Potential Visualizer indirectly through montage_at_points()."""
    f, _d, _t = _load_basis(sol_npz)
    currents_A = currents_mA_to_A(currents_mA, f.values.shape[0])
    return (currents_A[:, None] * f.values).sum(axis=0)


def montage_at_points(currents_mA, pts, sol_npz=None):
    """V (volts) at arbitrary (N,3) mm points in document coordinates; NaN
    at any point outside the solved FEM domain (TetField's own default)."""
    f, _d, _t = _load_basis(sol_npz)
    currents_A = currents_mA_to_A(currents_mA, f.values.shape[0])
    phi_at_pts = f(pts)  # (n_contacts, N)
    return (currents_A[:, None] * phi_at_pts).sum(axis=0)


def sample_tissue(tissue, currents_mA, sol_npz=None, surface_path=None):
    """Return (verts mm, tris, V volts) for one of config.TISSUE_STL's keys,
    sampled at the tissue's OWN STL vertices -- the real anatomy, not the
    synthetic gmsh box the FEM solve ran on.

    FOUND RUNNING THIS: TetField.locate()'s barycentric test (tolerance
    1e-9) fails for a real, non-negligible slice of these STL vertices --
    10-14% on dura/white/grey for the default bipolar case -- because the
    tet mesh's tissue-labelled boundary does not exactly coincide with the
    true STL surface (fem/README.md: "interfaces are resolved to element
    size, not followed exactly", and HANDOFF.md's own note that single-ray
    point classification can invert on a grazing edge). Measured: every
    such point sits within 0.5 mm of the nearest solved FEM node -- exactly
    the coarsest target element size in this mesh, never more -- so this is
    a resolution artifact at the boundary, not points genuinely outside the
    solved domain. Rather than leave gaps/garbage in the surface plot, these
    are filled from the nearest solved node's own value (nearest-neighbour,
    not extrapolated). Reported below, not hidden."""
    path = surface_path or C.TISSUE_BODIES[tissue]
    verts, tris = read_surface(path)
    V = montage_at_points(currents_mA, verts, sol_npz=sol_npz)
    nan_mask = np.isnan(V)
    n_nan = int(nan_mask.sum())
    if n_nan:
        _f, d, node_tree = _load_basis(sol_npz)
        V_nodes = montage_at_nodes(currents_mA, sol_npz=sol_npz)
        dist, idx = node_tree.query(verts[nan_mask])
        n_far = int((dist > FILL_TOL_MM).sum())
        V[nan_mask] = V_nodes[idx]
        print(
            "NOTE: %d/%d vertices of '%s' were outside TetField's strict "
            "containment test; filled from the nearest solved mesh node "
            "(max %.3f mm away, tolerance %.1f mm). %s"
            % (
                n_nan, len(V), tissue, float(dist.max()), FILL_TOL_MM,
                ("ALL within tolerance." if n_far == 0 else
                 "WARNING: %d of these were FURTHER than the fill tolerance "
                 "-- treat as genuinely suspect, not a resolution artifact."
                 % n_far),
            )
        )
    return verts, tris, V


def choose_clamp_v(currents_mA, sol_npz=None):
    """Pick a bulk-tissue-legible colour-legend clamp: the 95th percentile of
    |V| pooled over the three target tissues' own vertices, for THIS
    montage, rounded up to a tidy number a clinician can read at a glance."""
    pooled = []
    for t in CORE_TISSUES:
        _v, _t, V = sample_tissue(t, currents_mA, sol_npz=sol_npz)
        finite = V[np.isfinite(V)]
        if finite.size:
            pooled.append(np.abs(finite))
    if not pooled:
        return 100.0
    p95 = float(np.percentile(np.concatenate(pooled), 95))
    for tidy in (1, 2, 5, 10, 20, 25, 50, 75, 100, 150, 200, 250, 500, 750, 1000, 1500, 2000):
        if tidy >= p95:
            return float(tidy)
    return float(np.ceil(p95 / 500.0) * 500.0)


def write_vtp(path, verts, tris, V, clamp_v):
    """Write a vtkPolyData (.vtp) surface -- same triangulation as the source
    STL -- carrying both 'V_volts' (true) and 'V_clamped_<N>V' (saturated to
    +/-clamp_v) as point-data arrays. Loadable directly by FreeCAD's
    Fem::FemPostPipeline.read() (vtp is one of its supported extensions --
    checked against FemPostPipeline::canRead in this build's C++ source)."""
    import vtk
    from vtk.util import numpy_support as vnp

    verts = np.ascontiguousarray(verts, dtype=np.float64)
    tris = np.ascontiguousarray(tris, dtype=np.int64)

    pts = vtk.vtkPoints()
    pts.SetData(vnp.numpy_to_vtk(verts, deep=True))

    conn = np.empty((len(tris), 4), dtype=np.int64)
    conn[:, 0] = 3
    conn[:, 1:] = tris
    id_arr = vnp.numpy_to_vtkIdTypeArray(conn.ravel(), deep=True)
    polys = vtk.vtkCellArray()
    polys.SetCells(len(tris), id_arr)

    poly = vtk.vtkPolyData()
    poly.SetPoints(pts)
    poly.SetPolys(polys)

    v_raw = vnp.numpy_to_vtk(np.ascontiguousarray(V, dtype=np.float64), deep=True)
    v_raw.SetName("V_volts")
    poly.GetPointData().AddArray(v_raw)

    clamp_name = "V_clamped_%gV" % clamp_v
    v_clamp = np.clip(V, -clamp_v, clamp_v)
    v_clamp_arr = vnp.numpy_to_vtk(np.ascontiguousarray(v_clamp, dtype=np.float64), deep=True)
    v_clamp_arr.SetName(clamp_name)
    poly.GetPointData().AddArray(v_clamp_arr)

    poly.GetPointData().SetActiveScalars(clamp_name)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(path)
    writer.SetInputData(poly)
    writer.SetDataModeToAppended()
    writer.EncodeAppendedDataOff()
    writer.Write()
    return path


def surface_paths(sol_npz=None):
    """Use the selected run's immutable lead export, never the default lead."""
    paths = {t: tuple(C.TISSUE_BODIES[t]) for t in TISSUES}
    if sol_npz:
        lead_dir = os.path.join(os.path.dirname(sol_npz), 'lead')
        if os.path.isdir(lead_dir):
            contacts, insulator = C.lead_stls(lead_dir)
            paths.update({'contact_%d' % i: path for i, path in sorted(contacts.items())})
            import glob
            separate = sorted(glob.glob(os.path.join(lead_dir, 'insulator_lead_*.stl')))
            if separate:
                paths.update({'insulator_' + os.path.basename(p)[len('insulator_lead_'):-4]: p for p in separate})
            else:
                paths['insulator'] = insulator
    return paths


def run_montage(currents_mA, out_dir=None, clamp_v=None, sol_npz=None):
    """Write all anatomical and lead surfaces for the given per-contact currents
    (mA, + = source/anode, - = sink/cathode). Returns a stats dict for the
    caller (FreeCAD glue or the CLI below) to print/warn on."""
    out_dir = out_dir or os.path.join(C.OUT, "scs_viz")
    os.makedirs(out_dir, exist_ok=True)
    currents_mA = np.asarray(currents_mA, dtype=np.float64)
    if clamp_v is None:
        clamp_v = choose_clamp_v(currents_mA, sol_npz=sol_npz)
    stats = {
        "clamp_v": float(clamp_v),
        "charge_sum_mA": float(currents_mA.sum()),
        "tissues": {},
    }
    for t, surface_path in surface_paths(sol_npz).items():
        verts, tris, V = sample_tissue(t, currents_mA, sol_npz=sol_npz, surface_path=surface_path)
        path = os.path.join(out_dir, "%s.vtp" % t)
        write_vtp(path, verts, tris, V, clamp_v)
        finite = V[np.isfinite(V)]
        stats["tissues"][t] = dict(
            path=path,
            n_verts=int(len(V)),
            n_nan=int(np.isnan(V).sum()),
            v_min=float(finite.min()) if finite.size else float("nan"),
            v_max=float(finite.max()) if finite.size else float("nan"),
        )
    return stats


def self_check(sol_npz=None, atol=1e-3):
    """Reproduce fem/README.md's documented bipolar range from this module's
    OWN general montage() path -- not from a hand-copied constant -- and diff
    it node-for-node against solution.npz's own stored 'V' array (which is
    exactly phi[2]-phi[4], the field of record for that case). Returns
    (ok: bool, message: str).

    Only meaningful for a lead that actually HAS both config.SOURCE_CONTACT
    and config.SINK_CONTACT (true for the frozen 8-contact lead; also true
    for any live lead with >= SINK_CONTACT contacts, contiguously numbered
    from 1 -- see live_lead.find_live_lead()). For a lead with fewer
    contacts, assign_and_solve.solve_basis_fields() itself never computed a
    real bipolar 'V' for it (documented zero placeholder -- see that
    function's own docstring), so there is nothing to compare here: this
    returns ok=True with an explanatory message instead of a false FAIL,
    since a live N<5-contact lead is exactly the case this pipeline was
    just generalized to support, not an error."""
    f, d, _t = _load_basis(sol_npz)
    n_contacts = f.values.shape[0]
    if n_contacts < max(C.SOURCE_CONTACT, C.SINK_CONTACT):
        return True, (
            "self-check skipped: this lead has %d contact(s), fewer than "
            "config.SOURCE_CONTACT/SINK_CONTACT (%d/%d) -- no bipolar "
            "reference field applies to it (see solve_basis_fields()'s "
            "'V' placeholder note); its real per-contact fields are in "
            "'phi', not 'V'." % (n_contacts, C.SOURCE_CONTACT, C.SINK_CONTACT)
        )
    currents_mA = default_bipolar_mA(n_contacts)
    V_general = montage_at_nodes(currents_mA, sol_npz=sol_npz)
    V_stored = d["V"].astype(np.float64)
    max_diff = float(np.abs(V_general - V_stored).max())
    lo, hi = float(V_general.min()), float(V_general.max())
    msg = (
        "general montage vs solution.npz['V']: max|diff| = %.3e V; "
        "range %.2f to %.2f V (README documents -1369.91 to +1297.67 V "
        "for the frozen 8-contact lead)"
        % (max_diff, lo, hi)
    )
    return max_diff < atol, msg


def main():
    import argparse

    default_mA = default_bipolar_mA()
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        # nargs="+", not a fixed count: this CLI's --sol-npz can point at
        # any lead's solution.npz (a live N-contact one included), so the
        # required length is only known once that file is loaded --
        # currents_mA_to_A() validates it against the actual phi.shape[0]
        # below, with a clear error, rather than argparse enforcing a
        # number that's only even correct for the frozen 8-contact lead.
        "--currents-mA", type=float, nargs="+", default=None,
        help=(
            "Per-contact currents, mA (as many values as the target "
            "solution.npz has basis fields for -- 8 for the frozen "
            "fem/out/solution.npz), + = anode/source, - = cathode/sink; "
            "default = the validated bipolar case for that lead "
            "(contact %d = %+.0f mA, contact %d = %+.0f mA, for the frozen "
            "8-contact lead)"
            % (C.SOURCE_CONTACT, default_mA[C.SOURCE_CONTACT - 1],
               C.SINK_CONTACT, default_mA[C.SINK_CONTACT - 1])
        ),
    )
    p.add_argument("--clamp-v", type=float, default=None,
                    help="override the auto-chosen bulk-tissue clamp, volts")
    p.add_argument("--out-dir", default=None)
    p.add_argument("--sol-npz", default=None,
                    help="solution.npz to montage against; default = the "
                         "frozen fem/out/solution.npz")
    args = p.parse_args()

    ok, msg = self_check()
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        print(
            "REFUSING to proceed: general montage arithmetic disagrees with "
            "the field of record -- see fem/README.md RESULTS. Don't "
            "second-guess the physics, report this instead."
        )
        sys.exit(1)

    currents_mA = np.asarray(args.currents_mA) if args.currents_mA is not None else default_mA
    stats = run_montage(currents_mA, out_dir=args.out_dir, clamp_v=args.clamp_v)
    print(
        "clamp: +/-%.3g V (%s)"
        % (stats["clamp_v"], "user override" if args.clamp_v is not None else "auto, 95th pct |V|")
    )
    print("charge sum: %.6g mA" % stats["charge_sum_mA"])
    for t, s in stats["tissues"].items():
        print(
            "  %-6s %6d verts  V in [%8.2f, %8.2f] V  NaN=%d  -> %s"
            % (t, s["n_verts"], s["v_min"], s["v_max"], s["n_nan"], s["path"])
        )


if __name__ == "__main__":
    main()
