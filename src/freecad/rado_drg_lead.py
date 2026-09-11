#!/usr/bin/env python3
"""Place RADO's OWN DRG lead on any of the eight ganglia, by rigid transform.

WHAT CHANGED, AND WHY
---------------------
A DRG lead used to be generated the way a dorsal one is: a segment plan swept
along a measured corridor, with contact count, length, gap and diameter as
spinboxes. That is the right way to build an epidural lead, which is a straight
run down a channel that the model measures directly. It was the wrong way to
build a DRG lead, and the geometry showed it -- the sweep had to be nudged off
the corridor centre by a solved-for standoff to stop it cutting into the
ganglion, the insulator's curvature was whatever a polyline sweep gave, and the
result was a lead that fitted the foramen without ever being the lead RADO drew.

RADO's is the lead RADO drew. It ships in STL_files/ as five bodies -- four
platinum contacts and an insulator sheath -- hand-placed in SolidWorks against
the left third ganglion, and it is the only DRG lead in this model that anybody
has ever looked at and called correct. So this module stops generating DRG
geometry and starts COPYING it: the five meshes are transformed rigidly, vertex
by vertex, from the ganglion RADO built them against onto whichever ganglion is
asked for. Contact lengths, spacings, diameters and the sheath's curvature come
along unchanged, because a rigid transform cannot change them. That is the
point.

Dorsal and ventral leads are untouched. They are still swept and still fully
parametric; see make_scs_lead.py and build_lead_config.py.

WHICH GANGLION IS RADO'S OWN, AND HOW THAT WAS ESTABLISHED
-----------------------------------------------------------
L3 -- the patient's LEFT side, third ganglion counting rostral to caudal, which
in this model's own numbering (see lead_defaults.yaml's caveat: these are not
lumbar levels) sits at z = 93.4. Established by measurement, not assumption, and
the margin is not close:

    ganglion   |centroid - lead centroid|    nearest lead-surface-to-core gap
    L3                      3.28 mm                        0.29 mm
    L2                     28.16 mm                       23.22 mm
    R3                     33.48 mm                       22.10 mm
    ... the other five, 35 to 65 mm

An order of magnitude between first and second place on both measures, and the
whole lead lies at x = 67.1 .. 82.7, entirely on the +x (anatomical left) side
of the midline at x = 56.60. There is no ambiguity to resolve.

THE FRAME, AND WHY IT IS BUILT THE WAY IT IS
---------------------------------------------
Every ganglion gets the same three-step construction, so that "where RADO's lead
sits relative to L3" is a set of coordinates that can be re-read at any other
ganglion:

    origin  the ganglion core's AREA-WEIGHTED centroid. Area-weighted rather
            than a vertex mean because RADO tessellates curvature finely and a
            vertex mean drifts toward the curved end of a body.

    e1      the ganglion core's LONG AXIS -- the principal eigenvector of its
            area-weighted covariance -- oriented OUTWARD, away from the midline.
            This is "along the root": the DRG is a prolate football and its long
            axis is the direction the root runs.

    e2      anatomical DORSAL, Gram-Schmidt'd against e1 (global +Y minus its e1
            component, normalised).

    e3      anatomical ROSTRAL, Gram-Schmidt'd against e1 and e2.

WHY THE MESH AXIS AND NOT THE CORRIDOR TANGENT. measure_foramen.py fits a
quartic y(u), z(u) through each foramen and its derivative is an obvious
candidate for e1. It is not used, for a measured reason: the eight cores are
CONGRUENT CAD COPIES of two master bodies -- all four left cores are 5600
facets with identical sorted edge lengths to 3e-5 mm and signed volume
100.0718 mm3; all four right cores are 6020 facets and 101.0664 mm3 -- and
their principal axes come out as exact whole-degree rotations (0.978/-0.208 is
cos/sin 12 degrees, 0.966/-0.259 is 15 degrees). The mesh axis is recovering
the rotation RADO typed into SolidWorks. The quartic's derivative, fitted to
0.13-0.15 mm rms through ray-cast stations, disagrees with it by 1.8 to 9.7
degrees, and that disagreement is the fit's noise, not the anatomy's. The
corridor tangent is still COMPUTED and reported, as a cross-check that the two
descriptions of the foramen have not diverged; it just does not drive anything.

WHY e2 AND e3 COME FROM THE GLOBAL AXES AND NOT FROM THE MESH. The second and
third eigenvalues of every core are DEGENERATE -- 2.19 and 2.19 on the left,
2.20 and 2.20 on the right -- because the football is axisymmetric. Its
"second principal axis" is numerical noise and would spin the lead around the
root at random. So the roll is fixed externally, by the global dorsal and
rostral directions, which is both reproducible and anatomically what is meant:
a DRG lead goes into the superior aspect of the foramen whichever side it is on.

A useful consequence of the cores being congruent: for a SAME-SIDE move (L3 to
L1, say) the transform maps the source core exactly onto the target core, so the
lead's standoff from its ganglion is preserved to floating point. What differs
between targets is the foramen around it, and that is what the fit check
measures.

THE MIRROR, WHICH IS THE TRAP
------------------------------
Going from a left ganglion to a right one is a REFLECTION, not a rotation. The
naive fix -- rotate the lead 180 degrees about z to point it the other way --
also swaps dorsal for ventral, and puts the lead in front of the cord.

Here the reflection is not bolted on; it FALLS OUT of the frame, because the
frame is defined anatomically rather than by coordinates. On a left ganglion
e1 is roughly +x, e2 is +y and e3 works out as +z, and [e1 e2 e3] has
determinant +1. On a right ganglion e1 is roughly -x while e2 is still dorsal
and e3 is still rostral, and the determinant is -1. So

    M = R_target . R_source^T          det = +1 same side, -1 across sides

is a rotation within a side and an improper (reflection-bearing) transform
across sides, automatically and without a special case. It is checked, not
assumed: place() records the determinant and asserts it is +/-1 to 1e-12, and
for a cross-side move it factors M as (sagittal mirror) x (residual rotation)
and reports how far that residual is from the identity -- which is a direct
measurement of how asymmetric this model's left and right foramina are.

AND THE WINDING, WHICH IS THE TRAP INSIDE THE TRAP
---------------------------------------------------
A reflection reverses triangle orientation. Transform the three vertices of a
facet through a negative-determinant matrix and the facet's normal, recomputed
by the right-hand rule, now points INTO the solid. An STL like that is not
merely cosmetically wrong:

  - the point-in-mesh tests this pipeline is built on count ray crossings for
    parity (build_lead_config.inside), and inverted normals do not change
    crossing COUNTS -- so the containment check would still work, and would
    quietly stop being evidence of anything;
  - Simpleware and Ansys both read normals to decide which side of a surface is
    material, and an inside-out lead meshes as a void or fails outright.

So apply_transform() reverses the vertex order of every facet whenever the
determinant is negative, and the result is verified by signed volume: sum over
facets of v0 . (v1 x v2) / 6, which is positive for a closed mesh with outward
normals and negative for an inverted one. All five of RADO's STLs are watertight
(zero non-manifold edges) and positive, so the transformed copies must be too,
and with the same magnitude, because a rigid map preserves volume.

Note that FreeCAD's own Mesh.Mesh.Volume disagrees with this -- it reports
1.0729, 1.1979, 1.1823, 1.1562 for the four contacts, which are congruent
bodies with identical surface area and identical signed volume 1.09815. Whatever
it is doing, it is not a reliable orientation test here, so the check is done in
numpy where the arithmetic is visible.

WHAT IS STILL A PARAMETER
--------------------------
    target          which ganglion: L1-L4, R1-R4
    lateral_offset  slide along the root (+e1, further out), mm
    y_offset        nudge dorsally (+e2), mm
    z_offset        nudge rostrally (+e3), mm

All three offsets are in the TARGET'S OWN FRAME, which is what makes them mean
the same thing on both sides: +z_offset is rostral at every ganglion, left or
right, even though the transform that got there was a reflection.

What is NOT a parameter any more, for any DRG lead: contacts, contact_length,
gap, diameter, tail. They are RADO's, they are fixed, and build_lead_config
rejects them rather than ignoring them. They are still REPORTED -- measured back
off the placed geometry by hardware_spec(), in the lead's own frame -- so the
report can still say what the hardware is.

ganglion_clearance is also gone. It existed because a swept lead had no
inherent relationship to its ganglion and one had to be solved for; a copy of
RADO's lead has RADO's relationship built in. If a particular foramen needs the
lead held further off, z_offset is the knob, and the fit check prints the
standoff so you know how much to ask for.

NO FREECAD HERE
---------------
This module is numpy and struct. It reads and writes binary STL itself, so the
transform, the winding fix and the volume check can all be run and tested
without a FreeCAD interpreter:

    python3 src/freecad/rado_drg_lead.py --target R2 --show
    python3 src/freecad/rado_drg_lead.py --all

The containment and clearance checks in build_lead_config._validate_drg do need
FreeCAD, because they ray-cast the anatomy meshes.
"""

import argparse
import json
import os
import struct
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
DEFAULT_STL = os.path.join(REPO, "STL_files")
DEFAULT_FORAMEN = os.path.join(HERE, "foraminal_corridors.json")

MIDLINE_X = 56.60

# The ganglion RADO's own lead was built against; see the docstring for the
# measurement that establishes it. verify_native_target() re-runs that
# measurement from the meshes, so this constant can never quietly go stale.
NATIVE_TARGET = "L3"

# RADO's five hardware bodies, in the order they are exported. The FILENAMES
# ARE KEPT EXACTLY: ../ansys/tissue_map.yaml matches the object Names these
# sanitise to, so renaming a contact costs it its conductivity and its colour.
#
# The numbers are RADO's and are NOT in spatial order -- running out along the
# root the four contacts go 3, 1, 2, 4. They are deliberately not renumbered: a
# copy of a body should be that body. hardware_spec() reports the spatial order
# instead, and export writes it into the manifest, so nobody wires up a bipolar
# pair from the file numbers by mistake.
BODIES = (
    ("contact", 1, "SCS Lead Electrode 1.stl"),
    ("contact", 2, "SCS Lead Electrode 2.stl"),
    ("contact", 3, "SCS Lead Electrode 3.stl"),
    ("contact", 4, "SCS Lead Electrode 4.stl"),
    ("insulator", 1, "SCS Lead Insulator.stl"),
)

# Offsets, in the target ganglion's own frame. The only DRG parameters left.
OFFSET_PARAMS = ("lateral_offset", "y_offset", "z_offset")


# --------------------------------------------------------------------------
# binary STL, read and written here so this module needs no FreeCAD
# --------------------------------------------------------------------------
def read_stl(path):
    """An (N, 3, 3) float64 array of facet vertices. Binary or ASCII.

    The binary reader is find_floating_bodies.read_stl_vertices' -- the 50-byte
    record is 12 bytes of stored normal then 36 of vertices then 2 of attribute
    -- reshaped per facet rather than flat, because every operation here is
    per-facet. The STORED normal is deliberately discarded: it is redundant with
    the winding, the two disagree in plenty of real files, and this module's
    whole job is to keep the winding right and recompute normals from it.
    """
    with open(path, "rb") as fh:
        head = fh.read(84)
        if len(head) < 84:
            return np.empty((0, 3, 3), np.float64)
        ntri = struct.unpack("<I", head[80:84])[0]
        rest = fh.read()
    if ntri > 0 and len(rest) == ntri * 50:
        rec = np.frombuffer(rest, dtype=np.uint8).reshape(ntri, 50)
        v = rec[:, 12:48].copy().view("<f4").reshape(ntri, 3, 3)
        return np.asarray(v, dtype=np.float64)
    pts = []
    with open(path, "r", errors="ignore") as fh:
        for line in fh:
            if "vertex" in line:
                pts.append([float(v) for v in line.split()[1:4]])
    return np.asarray(pts, dtype=np.float64).reshape(-1, 3, 3)


def facet_normals(tris):
    """Unit facet normals FROM THE WINDING, right-hand rule. (N, 3).

    Degenerate facets (zero area) get a zero normal rather than a NaN, so one
    sliver in a source mesh cannot poison a whole file's normals.
    """
    n = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    ln = np.linalg.norm(n, axis=1)
    ok = ln > 0
    out = np.zeros_like(n)
    out[ok] = n[ok] / ln[ok, None]
    return out


def write_stl(path, tris, header=b"RADO-SCS DRG lead, rigidly transformed copy"):
    """Write a binary STL, with normals recomputed from the winding.

    Recomputed rather than carried through: after a reflection the source file's
    stored normals are wrong by construction, and a normal that contradicts its
    facet's winding is exactly the failure this module exists to prevent.
    """
    tris = np.ascontiguousarray(tris, dtype=np.float32)
    nrm = facet_normals(tris.astype(np.float64)).astype(np.float32)
    rec = np.zeros((len(tris), 50), np.uint8)
    rec[:, 0:12] = nrm.view(np.uint8).reshape(-1, 12)
    rec[:, 12:48] = tris.reshape(len(tris), 9).view(np.uint8).reshape(-1, 36)
    with open(path, "wb") as fh:
        fh.write(header[:80].ljust(80, b"\0"))
        fh.write(struct.pack("<I", len(tris)))
        fh.write(rec.tobytes())


def signed_volume(tris):
    """Divergence-theorem volume. Positive iff a closed mesh has OUTWARD normals.

    This is the orientation test. It does not care where the origin is (for a
    closed surface the origin cancels), it uses only the winding, and it changes
    sign -- not magnitude -- when the winding is reversed. See the docstring for
    why FreeCAD's own Mesh.Volume is not used for this.
    """
    return float(np.einsum("ij,ij->i", tris[:, 0], np.cross(tris[:, 1], tris[:, 2])).sum() / 6.0)


def surface_area(tris):
    return float((0.5 * np.linalg.norm(
        np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0]), axis=1)).sum())


def open_edges(tris, decimals=5):
    """How many edges are not shared by exactly two facets. 0 means watertight.

    Vertices are rounded before being matched, because an STL stores each
    facet's vertices independently and two facets meeting at a corner carry two
    separately-rounded float32 copies of it.
    """
    q = np.round(tris.reshape(-1, 3), decimals)
    _u, inv = np.unique(q, axis=0, return_inverse=True)
    inv = inv.reshape(-1, 3)
    e = np.concatenate([inv[:, [0, 1]], inv[:, [1, 2]], inv[:, [2, 0]]])
    _k, counts = np.unique(np.sort(e, axis=1), axis=0, return_counts=True)
    return int((counts != 2).sum())


def area_centroid(tris):
    """Area-weighted centroid, (3,). See measure_foramen.area_centroid."""
    a = 0.5 * np.linalg.norm(
        np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0]), axis=1)
    c = tris.mean(axis=1)
    return (c * a[:, None]).sum(0) / a.sum()


def principal_axes(tris):
    """Area-weighted principal axes of a mesh. (eigenvalues desc, columns).

    Weighted by facet area for the same reason the centroid is: an unweighted
    vertex covariance measures the tessellation as much as the shape.
    """
    a = 0.5 * np.linalg.norm(
        np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0]), axis=1)
    c = tris.mean(axis=1)
    o = (c * a[:, None]).sum(0) / a.sum()
    x = c - o
    cov = (x[:, :, None] * x[:, None, :] * a[:, None, None]).sum(0) / a.sum()
    w, v = np.linalg.eigh(cov)
    order = np.argsort(w)[::-1]
    return w[order], v[:, order]


# --------------------------------------------------------------------------
# the frame
# --------------------------------------------------------------------------
def _unit(v):
    n = float(np.linalg.norm(v))
    if n < 1e-12:
        raise ValueError("cannot normalise a zero-length direction")
    return np.asarray(v, float) / n


def polyval(coef, x):
    v = 0.0
    for c in coef:
        v = v * x + c
    return v


def polyder(coef):
    n = len(coef) - 1
    return [coef[i] * (n - i) for i in range(n)]


def load_foramen(path=DEFAULT_FORAMEN):
    if not os.path.exists(path):
        raise ValueError(
            "no foraminal corridor measurements at %s -- the ganglion frames are "
            "built from them. Run: freecadcmd src/freecad/measure_foramen.py"
            % os.path.relpath(path, REPO))
    with open(path) as fh:
        return json.load(fh)


def ganglion_frame(tag, foramen=None, stl_dir=DEFAULT_STL, cores=None):
    """The local frame at one ganglion. See the module docstring for the choices.

    Returns a dict with `origin` (3,), `basis` (3, 3) whose COLUMNS are
    e1/e2/e3, `det` (+1 left, -1 right), and the diagnostics that say the
    construction is well conditioned: the eigenvalue spread that justifies
    taking e1 from the mesh, the degeneracy that forbids taking e2 and e3 from
    it, and the angle to the corridor tangent as an independent cross-check.
    """
    if foramen is None:
        foramen = load_foramen()
    if tag not in foramen["ganglia"]:
        raise ValueError("no ganglion named %r in this model (have: %s)"
                         % (tag, ", ".join(sorted(foramen["ganglia"]))))
    info = foramen["ganglia"][tag]
    sign = int(info["x_sign"])

    if cores is None:
        cores = {}
    core = cores.get(info["drg_core_stl"])
    if core is None:
        core = read_stl(os.path.join(stl_dir, info["drg_core_stl"]))
        cores[info["drg_core_stl"]] = core

    origin = area_centroid(core)
    w, v = principal_axes(core)

    e1 = v[:, 0]
    if np.dot(e1, [sign, 0.0, 0.0]) < 0:        # orient OUTWARD, away from midline
        e1 = -e1
    e1 = _unit(e1)

    # e2 = dorsal, e3 = rostral, both Gram-Schmidt'd against what is already
    # fixed. The residuals are asserted to be substantial rather than assumed:
    # a near-zero one would mean the root runs along the axis being projected,
    # and the frame would be ill conditioned rather than wrong-looking.
    r2 = np.array([0.0, 1.0, 0.0]) - np.dot([0.0, 1.0, 0.0], e1) * e1
    if np.linalg.norm(r2) < 0.5:
        raise ValueError("%s: the root axis is within 60 degrees of dorsal; the "
                         "frame's roll cannot be fixed from the global axes" % tag)
    e2 = _unit(r2)
    r3 = np.array([0.0, 0.0, 1.0]) - np.dot([0.0, 0.0, 1.0], e1) * e1 \
        - np.dot([0.0, 0.0, 1.0], e2) * e2
    if np.linalg.norm(r3) < 0.5:
        raise ValueError("%s: rostral lies in the plane of the root axis and "
                         "dorsal; the frame is degenerate" % tag)
    e3 = _unit(r3)

    basis = np.column_stack([e1, e2, e3])
    det = float(np.linalg.det(basis))
    # Orthonormality is a property of the construction, so a failure here is a
    # bug in it rather than a property of the anatomy. Checked anyway: every
    # claim this module makes about reflections rests on it.
    err = float(np.abs(basis.T @ basis - np.eye(3)).max())
    if err > 1e-9:
        raise AssertionError("%s: frame is not orthonormal (max error %.2e)" % (tag, err))
    if abs(abs(det) - 1.0) > 1e-9:
        raise AssertionError("%s: frame determinant %.12f is not +/-1" % (tag, det))
    expect = 1.0 if info["side"] == "left" else -1.0
    if det * expect < 0:
        raise AssertionError(
            "%s is on the %s and its frame determinant is %+.1f; the whole "
            "left/right mirror rests on that being %+.1f"
            % (tag, info["side"], det, expect))

    # The independent description of the same direction: the tangent of the
    # fitted foraminal corridor at the ganglion's centre. Not used to build
    # anything -- reported so the two cannot silently diverge.
    u = float(info["ganglion_u_centre"])
    tangent = _unit([sign,
                     polyval(polyder(info["y_poly_coef_high_to_low"]), u),
                     polyval(polyder(info["z_poly_coef_high_to_low"]), u)])
    angle = float(np.degrees(np.arccos(np.clip(np.dot(e1, tangent), -1.0, 1.0))))

    return {
        "tag": tag,
        "side": info["side"],
        "x_sign": sign,
        "origin": origin,
        "basis": basis,
        "det": det,
        "core_stl": info["drg_core_stl"],
        "core_facets": len(core),
        "core_signed_volume": signed_volume(core),
        "eigenvalues": [float(x) for x in w],
        "axis_degeneracy": float(abs(w[1] - w[2]) / w[0]),
        "corridor_tangent": tangent,
        "axis_vs_corridor_deg": angle,
        "ganglion_u_centre": u,
        "usable_u": info["usable_u"],
    }


def frame_matrix(frame):
    """The 4x4 that takes frame coordinates to model coordinates."""
    m = np.eye(4)
    m[:3, :3] = frame["basis"]
    m[:3, 3] = frame["origin"]
    return m


def transform_between(src, dst, offsets=None):
    """The rigid map from the source ganglion's frame to the target's. 4x4.

    p -> dst.origin + M (p - src.origin) + (offsets in the TARGET frame),
    with M = R_dst R_src^T. det(M) is +1 within a side and -1 across sides,
    which is the reflection, and it is asserted rather than hoped for.

    The offsets are folded into the translation, so what comes out is still one
    rigid transform and not a transform plus a fudge applied afterwards.
    """
    if src["tag"] == dst["tag"]:
        # R R^T IS the identity -- exactly, by definition, not to within
        # rounding. Computing the product would put 1e-16 off the diagonal and
        # move every vertex of a lead that has not been asked to move, which
        # would mean "copy RADO's lead onto its own ganglion" did not return
        # RADO's lead. It is the one assertion this whole mechanism can make
        # bit-exactly, so it is not thrown away for the sake of uniformity.
        m = np.eye(3)
    else:
        m = dst["basis"] @ src["basis"].T
    det = float(np.linalg.det(m))
    if abs(abs(det) - 1.0) > 1e-9:
        raise AssertionError("the ganglion-to-ganglion map has determinant "
                             "%.12f; it is not rigid" % det)
    shift = np.zeros(3)
    if offsets:
        shift = (float(offsets.get("lateral_offset", 0.0)) * dst["basis"][:, 0]
                 + float(offsets.get("y_offset", 0.0)) * dst["basis"][:, 1]
                 + float(offsets.get("z_offset", 0.0)) * dst["basis"][:, 2])
    t = np.eye(4)
    t[:3, :3] = m
    t[:3, 3] = dst["origin"] + shift - m @ src["origin"]
    return t


def apply_transform(tris, matrix):
    """Transform facet vertices, REVERSING THE WINDING when the map reflects.

    The whole reason this function exists rather than a one-line matrix
    multiply. A negative-determinant map sends a right-handed triple to a
    left-handed one, so every facet's right-hand-rule normal ends up pointing
    into the solid. Reversing the vertex order of each facet -- (v0, v1, v2)
    becomes (v2, v1, v0) -- restores it. Verified by signed_volume(), which the
    callers check rather than trust.
    """
    r = np.asarray(matrix)[:3, :3]
    t = np.asarray(matrix)[:3, 3]
    out = tris @ r.T + t
    if np.linalg.det(r) < 0:
        out = out[:, ::-1, :]
    return np.ascontiguousarray(out)


def sagittal_mirror():
    """The exact reflection in the sagittal plane x = MIDLINE_X, as a 4x4."""
    m = np.eye(4)
    m[0, 0] = -1.0
    m[0, 3] = 2.0 * MIDLINE_X
    return m


def sagittal_decomposition(matrix, src, dst):
    """Factor a cross-side map as (sagittal mirror) then (rigid correction).

    A cross-side placement is a reflection plus whatever it takes to land on the
    target, and separating the two says how asymmetric this model actually is:
    if the left and right foramina were perfect mirror images, mirroring RADO's
    lead about x = MIDLINE_X would put it on the target with NOTHING left to
    correct.

    Two numbers describe what is left. The residual ROTATION is read off the
    correction directly. The residual DISPLACEMENT is deliberately NOT the
    correction's translation column -- that is measured about the model origin,
    40-odd millimetres away, so it is dominated by the rotation's lever arm and
    means nothing anatomically. It is measured where it matters instead: how far
    the mirror image of the source ganglion's centre lands from the target
    ganglion's centre.
    """
    mirror = sagittal_mirror()
    resid = np.linalg.inv(mirror) @ np.asarray(matrix)
    r = resid[:3, :3]
    ang = float(np.degrees(np.arccos(np.clip((np.trace(r) - 1.0) / 2.0, -1.0, 1.0))))
    mirrored_origin = (mirror @ np.append(src["origin"], 1.0))[:3]
    off = dst["origin"] - mirrored_origin
    return {
        "residual_rotation_deg": ang,
        "residual_det": float(np.linalg.det(r)),
        "ganglion_asymmetry_mm": float(np.linalg.norm(off)),
        "ganglion_asymmetry": off,
    }


# --------------------------------------------------------------------------
# the lead
# --------------------------------------------------------------------------
def load_rado_bodies(stl_dir=DEFAULT_STL, cache=None):
    """RADO's five hardware meshes as [(kind, idx, filename, tris)].

    Read once and cached by the caller, because a set of leads reads the same
    five files for every one of them.
    """
    if cache is not None and "rado_bodies" in cache:
        return cache["rado_bodies"]
    out = []
    for kind, idx, fn in BODIES:
        path = os.path.join(stl_dir, fn)
        if not os.path.exists(path):
            raise ValueError(
                "RADO's DRG lead is missing from %s: %s. A DRG lead is a copy of "
                "those five files, so without them there is nothing to place."
                % (os.path.relpath(stl_dir, REPO), fn))
        out.append((kind, idx, fn, read_stl(path)))
    if cache is not None:
        cache["rado_bodies"] = out
    return out


def _body_axis(tris, prefer):
    """A single body's own long axis, oriented to agree with `prefer`."""
    _w, v = principal_axes(tris)
    e = v[:, 0]
    return -e if np.dot(e, prefer) < 0 else e


def _cylinder_metrics(tris, axis):
    """(length, diameter) of a body about `axis`, measured from the body itself.

    Length is the extent along the axis; diameter is twice the largest distance
    from the axis line through the body's own area centroid. Both are measured
    about the BODY's axis rather than the lead's or the ganglion frame's,
    because a contact is tilted a few degrees relative to both and projecting a
    tilted cylinder onto someone else's axis reads long in length and wildly
    long in diameter -- an early version of this function reported RADO's 1.25
    mm contacts as 7.8 mm across for exactly that reason.
    """
    v = tris.reshape(-1, 3) - area_centroid(tris)
    s = v @ axis
    radial = np.linalg.norm(v - np.outer(s, axis), axis=1)
    return float(s.max() - s.min()), 2.0 * float(radial.max())


def _banded_diameter(tris, axis, bands=40):
    """Diameter of a CURVED body, measured band by band and taken as the median.

    A single radius about one straight axis is meaningless for the sheath, which
    follows the foramen's curve: the far end of the curve is several millimetres
    off the axis and would be reported as the diameter. Slicing along the axis
    and measuring each thin band about its OWN centre measures the cross-section
    instead. The median across bands ignores the rounded end caps, where a band
    is genuinely narrower.
    """
    v = tris.reshape(-1, 3)
    s = v @ axis
    edges = np.linspace(s.min(), s.max(), bands + 1)
    out = []
    for i in range(bands):
        sel = v[(s >= edges[i]) & (s < edges[i + 1])]
        if len(sel) < 8:
            continue
        p = sel - sel.mean(axis=0)
        p = p - np.outer(p @ axis, axis)
        out.append(2.0 * float(np.linalg.norm(p, axis=1).max()))
    return float(np.median(out)) if out else 0.0


def hardware_spec(bodies):
    """Measure the hardware back off the geometry. Nothing here is asserted.

    Contact count, contact length, gap, pitch, diameter and tail are no longer
    parameters -- they are RADO's, fixed by the meshes. They are still worth
    SAYING, so they are measured, and measured in the way each quantity is
    actually defined:

        pitch          centre-to-centre distance between consecutive contacts.
                       The most robust number here, because it depends on no
                       axis at all.
        contact length and diameter
                       about EACH CONTACT'S OWN principal axis, which is the
                       only axis a cylinder's length and diameter mean anything
                       about.
        gap            pitch minus length, i.e. the insulator showing between
                       two neighbours.
        sheath diameter
                       banded, because the sheath is curved (see
                       _banded_diameter).
        tail           insulator beyond the end contacts, SIGNED. RADO's comes
                       out negative at both ends: the end contacts overhang the
                       sheath rather than sitting inboard of it, which is why
                       the parameters that used to reproduce this lead passed
                       --tail 0.

    `spatial_order` is RADO's electrode numbers sorted medial to lateral. They
    are not in file order -- they run 3, 1, 2, 4 -- and a bipolar pair chosen
    from the file numbers would not be the pair the clinician meant.

    NO FRAME IS NEEDED, deliberately. Every quantity here is intrinsic to the
    five meshes, so this can be called on RADO's unplaced originals -- which is
    what resolve_lead() does, to fill in the contact count and dimensions a DRG
    lead no longer takes as parameters, without having to read
    foraminal_corridors.json first. "Outward" is settled by the lead's own
    position relative to the midline rather than by a ganglion's axis.
    """
    contacts = [(kind, idx, fn, tris) for kind, idx, fn, tris in bodies
                if kind == "contact"]
    ins = next((t for k, _i, _f, t in bodies if k == "insulator"), None)

    if not contacts:
        return {"contacts": 0}
    cen = {idx: area_centroid(tris) for _k, idx, _f, tris in contacts}
    # Outward = away from the sagittal plane, in the direction this lead already
    # lies. True on either side and for a lead that has not been placed yet.
    all_x = np.concatenate([t.reshape(-1, 3)[:, 0] for _k, _i, _f, t in bodies])
    e1 = np.array([1.0 if all_x.mean() >= MIDLINE_X else -1.0, 0.0, 0.0])
    order = sorted(cen, key=lambda i: float(cen[i] @ e1))

    # The lead's own axis: the total-least-squares line through the contact
    # centroids, which describes where the array runs better than either the
    # ganglion's long axis or any single contact's. SVD rather than the
    # area-weighted covariance used for a mesh -- these are four bare points,
    # they carry no area, and feeding them to principal_axes() as degenerate
    # facets divides by a zero total area.
    pts = np.array([cen[i] for i in order])
    if len(pts) > 2:
        lead_axis = _unit(np.linalg.svd(pts - pts.mean(axis=0))[2][0])
    else:
        lead_axis = _unit(pts[-1] - pts[0])
    if np.dot(lead_axis, e1) < 0:
        lead_axis = -lead_axis

    lens, diams = [], []
    for _k, idx, _f, tris in contacts:
        ax = _body_axis(tris, lead_axis)
        ln, dm = _cylinder_metrics(tris, ax)
        lens.append(ln)
        diams.append(dm)
    pitches = [float(np.linalg.norm(pts[i + 1] - pts[i])) for i in range(len(pts) - 1)]

    spec = {
        "contacts": len(contacts),
        "contact_length": float(np.mean(lens)),
        "contact_length_range": [float(min(lens)), float(max(lens))],
        "diameter": float(np.mean(diams)),
        "diameter_range": [float(min(diams)), float(max(diams))],
        "pitch": float(np.mean(pitches)) if pitches else 0.0,
        "pitch_range": [float(min(pitches)), float(max(pitches))] if pitches else [0.0, 0.0],
        "gap": float(np.mean(pitches) - np.mean(lens)) if pitches else 0.0,
        "spatial_order": list(order),
        "array_mm": float(np.linalg.norm(pts[-1] - pts[0]) + np.mean(lens)),
        "array_centre": (pts[0] + pts[-1]) * 0.5,
    }
    if ins is not None:
        spec["insulator_diameter"] = _banded_diameter(ins, lead_axis)
        s_ins = (ins.reshape(-1, 3)) @ lead_axis
        c_lo = min(float((t.reshape(-1, 3) @ lead_axis).min()) for _k, _i, _f, t in contacts)
        c_hi = max(float((t.reshape(-1, 3) @ lead_axis).max()) for _k, _i, _f, t in contacts)
        spec["tail_medial"] = float(c_lo - s_ins.min())
        spec["tail_lateral"] = float(s_ins.max() - c_hi)
        spec["tail"] = float(min(spec["tail_medial"], spec["tail_lateral"]))
        spec["total_mm"] = float(max(s_ins.max(), c_hi) - min(s_ins.min(), c_lo))
    return spec


def place(target, foramen=None, stl_dir=DEFAULT_STL, lateral_offset=0.0,
          y_offset=0.0, z_offset=0.0, cache=None):
    """Copy RADO's DRG lead onto `target`. The one function this module is for.

    Returns a dict carrying the transformed meshes, the transform that produced
    them, both frames, the measured hardware, and -- the part that matters --
    the ORIENTATION AUDIT: per body, the signed volume before and after, whether
    the winding was reversed, and whether the volume was preserved in magnitude
    and sign. A caller that exports without looking at `orientation_ok` is
    shipping meshes it has not checked.
    """
    if cache is None:
        cache = {}
    if foramen is None:
        foramen = load_foramen()
    cores = cache.setdefault("cores", {})
    src = ganglion_frame(NATIVE_TARGET, foramen, stl_dir, cores)
    dst = ganglion_frame(target, foramen, stl_dir, cores)
    offsets = {"lateral_offset": lateral_offset, "y_offset": y_offset,
               "z_offset": z_offset}
    matrix = transform_between(src, dst, offsets)
    det = float(np.linalg.det(matrix[:3, :3]))
    mirrored = det < 0

    bodies, audit, ok = [], [], True
    for kind, idx, fn, tris in load_rado_bodies(stl_dir, cache):
        moved = apply_transform(tris, matrix)
        v0, v1 = signed_volume(tris), signed_volume(moved)
        a0, a1 = surface_area(tris), surface_area(moved)
        # Three independent things, all of which must hold:
        #   the winding still encloses a positive volume (outward normals);
        #   the magnitude is preserved (the map really was rigid);
        #   the area is preserved (ditto, and catches a scale slipping in).
        row = {
            "body": fn,
            "kind": kind,
            "index": idx,
            "facets": len(tris),
            "signed_volume_before": v0,
            "signed_volume_after": v1,
            "area_before": a0,
            "area_after": a1,
            "winding_reversed": bool(mirrored),
            "normals_outward": v1 > 0,
            "volume_preserved": abs(abs(v1) - abs(v0)) <= 1e-6 * max(1.0, abs(v0)),
            "area_preserved": abs(a1 - a0) <= 1e-6 * max(1.0, a0),
            "open_edges": open_edges(moved),
        }
        row["ok"] = bool(row["normals_outward"] and row["volume_preserved"]
                         and row["area_preserved"] and row["open_edges"] == 0)
        ok = ok and row["ok"]
        audit.append(row)
        bodies.append((kind, idx, fn, moved))

    out = {
        "target": target,
        "source": NATIVE_TARGET,
        "side": dst["side"],
        "cross_side": src["side"] != dst["side"],
        "mirrored": mirrored,
        "det": det,
        "matrix": matrix,
        "source_frame": src,
        "target_frame": dst,
        "offsets": offsets,
        "bodies": bodies,
        "orientation": audit,
        "orientation_ok": ok,
        "hardware": hardware_spec(bodies),
    }
    if mirrored != out["cross_side"]:
        raise AssertionError(
            "%s -> %s is %s a cross-side move but the transform %s a reflection "
            "(det %+.6f). The frame convention has broken."
            % (NATIVE_TARGET, target, "" if out["cross_side"] else "not",
               "is" if mirrored else "is not", det))
    if mirrored:
        out["sagittal"] = sagittal_decomposition(matrix, src, dst)
    return out


def vertices(placed, stride=1):
    """Every vertex of a placed lead, as an (N, 3) array. Deterministic.

    Sub-sampled by a fixed stride rather than randomly, so the same placement
    always reports the same numbers -- the same rule sample_points() follows.
    """
    v = np.vstack([b[3].reshape(-1, 3) for b in placed["bodies"]])
    return v[::stride] if stride > 1 else v


def sample_vertices(placed, max_samples=6000):
    v = vertices(placed)
    if len(v) <= max_samples:
        return v
    return v[:: len(v) // max_samples + 1]


# --------------------------------------------------------------------------
# exact point-to-surface distance, for the standoff
# --------------------------------------------------------------------------
def min_distance_to_surface(points, tris, chunk=256):
    """Smallest distance from any of `points` to the SURFACE of `tris`. mm.

    Point-to-triangle, not point-to-vertex. The distinction is the whole reason
    this exists: build_lead_config.min_distance_to_mesh compares vertices to
    vertices and therefore reads LONG -- by up to half a facet's width -- which
    is fine for the gauge it is used as and not fine for a number that is being
    quoted as "the standoff between the lead and the ganglion" and compared
    between targets.

    Triangles are pre-filtered against a vertex-to-vertex upper bound first,
    because the ganglion is 5600 facets and the lead 35000 vertices and the full
    product is not worth computing to find a gap that involves a few dozen of
    each.
    """
    points = np.asarray(points, float)
    tris = np.asarray(tris, float)
    if len(points) == 0 or len(tris) == 0:
        return float("inf")

    # Upper bound from a coarse vertex-to-vertex pass, then keep only the
    # triangles whose bounding box could possibly beat it.
    pv = points[:: max(1, len(points) // 400)]
    tv = tris.reshape(-1, 3)[:: max(1, len(tris) * 3 // 400)]
    bound = float(np.sqrt(((pv[:, None, :] - tv[None, :, :]) ** 2).sum(-1)).min())
    lo, hi = tris.min(axis=1), tris.max(axis=1)
    plo, phi = points.min(axis=0), points.max(axis=0)
    d = np.maximum(np.maximum(plo - hi, lo - phi), 0.0)
    near = np.linalg.norm(d, axis=1) <= bound + 1e-6
    tris = tris[near]
    if len(tris) == 0:
        return bound

    a, b, c = tris[:, 0], tris[:, 1], tris[:, 2]
    ab, ac = b - a, c - a
    best = bound
    for i in range(0, len(points), chunk):
        p = points[i:i + chunk][:, None, :]
        ap = p - a[None]
        d1 = (ab[None] * ap).sum(-1)
        d2 = (ac[None] * ap).sum(-1)
        bp = p - b[None]
        d3 = (ab[None] * bp).sum(-1)
        d4 = (ac[None] * bp).sum(-1)
        cp = p - c[None]
        d5 = (ab[None] * cp).sum(-1)
        d6 = (ac[None] * cp).sum(-1)
        va = d3 * d6 - d5 * d4
        vb = d5 * d2 - d1 * d6
        vc = d1 * d4 - d3 * d2

        denom = va + vb + vc
        denom = np.where(np.abs(denom) < 1e-30, 1e-30, denom)
        v = vb / denom
        w = vc / denom
        # Ericson's cascade, as a select. The interior case is the default and
        # each earlier condition overrides it, in the order the regions are
        # disjoint in.
        conds = [
            (d1 <= 0) & (d2 <= 0),                                   # vertex a
            (d3 >= 0) & (d4 <= d3),                                  # vertex b
            (d6 >= 0) & (d5 <= d6),                                  # vertex c
            (vc <= 0) & (d1 >= 0) & (d3 <= 0),                       # edge ab
            (vb <= 0) & (d2 >= 0) & (d6 <= 0),                       # edge ac
            (va <= 0) & ((d4 - d3) >= 0) & ((d5 - d6) >= 0),         # edge bc
        ]
        zero = np.zeros_like(d1)
        one = np.ones_like(d1)
        t_ab = np.where(np.abs(d1 - d3) < 1e-30, zero, d1 / np.where(np.abs(d1 - d3) < 1e-30, one, d1 - d3))
        t_ac = np.where(np.abs(d2 - d6) < 1e-30, zero, d2 / np.where(np.abs(d2 - d6) < 1e-30, one, d2 - d6))
        bcd = (d4 - d3) + (d5 - d6)
        t_bc = np.where(np.abs(bcd) < 1e-30, zero, (d4 - d3) / np.where(np.abs(bcd) < 1e-30, one, bcd))
        vs = [zero, one, zero, t_ab, zero, one - t_bc]
        ws = [zero, zero, one, zero, t_ac, t_bc]
        vv = np.select(conds, vs, default=v)
        ww = np.select(conds, ws, default=w)
        q = a[None] + vv[..., None] * ab[None] + ww[..., None] * ac[None]
        dist = np.linalg.norm(p - q, axis=-1).min()
        if dist < best:
            best = float(dist)
    return best


# --------------------------------------------------------------------------
# the check that NATIVE_TARGET is not stale
# --------------------------------------------------------------------------
def verify_native_target(foramen=None, stl_dir=DEFAULT_STL):
    """Re-derive which ganglion RADO's lead sits on. Returns (tag, ranking).

    The docstring's table, computed rather than quoted. Nearest ganglion by
    centroid distance, with the true surface standoff for the top few, so the
    margin between first and second place is visible and NATIVE_TARGET cannot
    quietly go stale if the meshes ever change.
    """
    if foramen is None:
        foramen = load_foramen()
    bodies = load_rado_bodies(stl_dir)
    lead = np.vstack([b[3].reshape(-1, 3) for b in bodies])
    centre = lead.mean(axis=0)
    rows = []
    for tag in sorted(foramen["ganglia"]):
        info = foramen["ganglia"][tag]
        core = read_stl(os.path.join(stl_dir, info["drg_core_stl"]))
        d = float(np.linalg.norm(area_centroid(core) - centre))
        rows.append({"tag": tag, "side": info["side"], "centroid_gap_mm": d,
                     "core": core})
    rows.sort(key=lambda r: r["centroid_gap_mm"])
    for r in rows[:3]:
        r["surface_gap_mm"] = min_distance_to_surface(lead[::5], r.pop("core"))
    for r in rows:
        r.pop("core", None)
    return rows[0]["tag"], rows


# --------------------------------------------------------------------------
# CLI -- reporting only; building a lead into a document is build_lead_config's
# --------------------------------------------------------------------------
def describe_placement(placed):
    """A placement as lines of text, for a report or the console."""
    f = placed["target_frame"]
    s = placed["source_frame"]
    hw = placed["hardware"]
    lines = [
        "RADO's DRG lead, copied from %s onto %s (%s side)"
        % (placed["source"], placed["target"], placed["side"]),
        "",
        "   source frame %-4s origin (%7.2f,%7.2f,%7.2f)  e1 (%+.4f,%+.4f,%+.4f)  det %+.1f"
        % (s["tag"], s["origin"][0], s["origin"][1], s["origin"][2],
           s["basis"][0, 0], s["basis"][1, 0], s["basis"][2, 0], s["det"]),
        "   target frame %-4s origin (%7.2f,%7.2f,%7.2f)  e1 (%+.4f,%+.4f,%+.4f)  det %+.1f"
        % (f["tag"], f["origin"][0], f["origin"][1], f["origin"][2],
           f["basis"][0, 0], f["basis"][1, 0], f["basis"][2, 0], f["det"]),
        "   root axis vs fitted corridor tangent: %.2f deg at %s, %.2f deg at %s"
        % (s["axis_vs_corridor_deg"], s["tag"], f["axis_vs_corridor_deg"], f["tag"]),
        "   transform determinant %+.9f -- %s"
        % (placed["det"],
           "REFLECTION (cross-side), windings reversed" if placed["mirrored"]
           else "rotation (same side), windings unchanged"),
    ]
    if placed["mirrored"]:
        sg = placed["sagittal"]
        lines.append(
            "   factored as the sagittal mirror at x = %.2f followed by a "
            "%.2f deg rotation" % (MIDLINE_X, sg["residual_rotation_deg"]))
        lines.append(
            ("   %s and %s are %.2f mm from being exact mirror images of each other"
             if s["tag"][1:] == f["tag"][1:] else
             "   the mirror image of %s's centre lands %.2f mm from %s's -- mostly "
             "the change of level, since these are not counterpart ganglia")
            % ((s["tag"], f["tag"], sg["ganglion_asymmetry_mm"])
               if s["tag"][1:] == f["tag"][1:]
               else (s["tag"], sg["ganglion_asymmetry_mm"], f["tag"])))
    off = placed["offsets"]
    if any(abs(v) > 1e-9 for v in off.values()):
        lines.append("   offsets in the target frame: %.2f mm along the root, "
                     "%.2f mm dorsal, %.2f mm rostral"
                     % (off["lateral_offset"], off["y_offset"], off["z_offset"]))
    lines += [
        "",
        "   hardware, measured back off the placed geometry (not parameters):",
        "      %d contacts, %.3f mm long, %.3f mm apart (%.3f mm pitch), %.3f mm across"
        % (hw["contacts"], hw["contact_length"], hw["gap"], hw["pitch"], hw["diameter"]),
        "      array %.2f mm, whole lead %.2f mm, sheath %.3f mm across, "
        "tails %+.2f / %+.2f mm"
        % (hw["array_mm"], hw.get("total_mm", float("nan")),
           hw.get("insulator_diameter", float("nan")),
           hw.get("tail_medial", float("nan")), hw.get("tail_lateral", float("nan"))),
        "      RADO's electrode numbers run %s medially to laterally"
        % " -> ".join(str(i) for i in hw["spatial_order"]),
        "",
        "   orientation audit (a reflection inverts normals unless the winding "
        "is reversed):",
        "      %-26s %10s %10s %9s %8s %s"
        % ("body", "vol before", "vol after", "reversed", "outward", "open edges"),
    ]
    for r in placed["orientation"]:
        lines.append("      %-26s %+10.5f %+10.5f %9s %8s %10d   %s"
                     % (r["body"], r["signed_volume_before"], r["signed_volume_after"],
                        "yes" if r["winding_reversed"] else "no",
                        "yes" if r["normals_outward"] else "NO",
                        r["open_edges"], "OK" if r["ok"] else "FAILED"))
    lines.append("      verdict: %s"
                 % ("all five bodies watertight, rigid and outward-facing"
                    if placed["orientation_ok"] else "FAILED -- do not export"))
    return lines


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--target", default=NATIVE_TARGET)
    ap.add_argument("--all", action="store_true", help="every ganglion")
    ap.add_argument("--native", action="store_true",
                    help="re-derive which ganglion RADO's own lead sits on")
    ap.add_argument("--lateral-offset", type=float, default=0.0)
    ap.add_argument("--y-offset", type=float, default=0.0)
    ap.add_argument("--z-offset", type=float, default=0.0)
    ap.add_argument("--stl-dir", default=DEFAULT_STL)
    ap.add_argument("--show", action="store_true")
    args = ap.parse_args()

    foramen = load_foramen()
    out = []
    if args.native:
        tag, rows = verify_native_target(foramen, args.stl_dir)
        out += ["RADO's own DRG lead sits on ganglion %s." % tag, "",
                "   %-4s %-6s %18s %18s" % ("tag", "side", "centroid gap", "surface gap")]
        for r in rows:
            out.append("   %-4s %-6s %15.2f mm %15s"
                       % (r["tag"], r["side"], r["centroid_gap_mm"],
                          ("%.2f mm" % r["surface_gap_mm"]) if "surface_gap_mm" in r else "-"))
        out.append("")

    targets = sorted(foramen["ganglia"]) if args.all else [args.target]
    cache = {}
    for tag in targets:
        placed = place(tag, foramen, args.stl_dir, args.lateral_offset,
                       args.y_offset, args.z_offset, cache)
        out += describe_placement(placed) + [""]
    sys.stdout.write("\n".join(out) + "\n")
    return 0


def run_as_freecadcmd_script():
    """True when a FreeCAD interpreter was handed THIS file; see apply_labels.py."""
    if len(sys.argv) < 2 or not os.path.basename(sys.argv[0]).lower().startswith("freecad"):
        return False
    return os.path.abspath(sys.argv[1]) == os.path.abspath(__file__)


if __name__ == "__main__":
    sys.exit(main())
elif run_as_freecadcmd_script():
    main()
