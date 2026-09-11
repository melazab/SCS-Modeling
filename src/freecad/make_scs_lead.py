"""Generate a parametric percutaneous SCS lead that follows the epidural canal.

Why this exists
---------------
The lead RADO ships is not the one we want to study. Measured from the geometry
(src/freecad/check_laterality.py) it sits 11-25 mm LEFT of midline at z = 93-97,
draped over the left DRG column -- a DRG lead -- on a curved insulator with only
4 contacts. Khadka et al. Figure 1 shows the clinical arrangement instead: a
percutaneous lead with 8 cylindrical contacts in the epidural space at midline.

This builds that lead at any rostro-caudal level, on either the dorsal or the
ventral side, which is the sweep the dorsal-vs-ventral study needs.

Following the canal matters
---------------------------
The T8-T10 spine here is kyphotic, so the epidural space migrates posteriorly
with height: at midline the dorsal channel centre runs y ~= 81.3 at z = 95 to
y ~= 88.4 at z = 145, about 7 mm over the length of an 8-contact lead. A
straight lead at fixed y would walk out of the space and through the dura. So
the lead is swept along the measured centreline y(z) from
src/freecad/epidural_corridor.json (produced by src/freecad/measure_corridor.py, which ray
casts the epidural mesh; cubic fit, 0.018 mm RMS over 98 mm).

Measured channel thickness at midline: dorsal 2.27-2.35 mm, ventral 1.67-1.74
mm. A 1.3 mm clinical lead fits either side; the script warns if the requested
diameter does not.

Frame (src/freecad/check_laterality.py): +X anatomical LEFT, midline x = 56.60 mm;
+Y posterior/DORSAL; +Z rostral. Millimetres throughout.

Construction: the lead is split along its length into alternating insulator and
contact segments, each a cylinder placed along the local tangent of the
centreline, so neighbours share a flat interface and never overlap -- which is
what the mesher and bonded electric contact want, and matches how RADO models
its own contacts (short full-diameter segments).

Usage -- run headless with FreeCAD:

    freecadcmd src/freecad/make_scs_lead.py
    MAKE_LEAD_ARGS="--side ventral --z-center 110" freecadcmd src/freecad/make_scs_lead.py
    MAKE_LEAD_ARGS="--contacts 8 --z-center 120 --tag sweep_z120" freecadcmd src/freecad/make_scs_lead.py

NOTE on freecadcmd: it *imports* this file rather than running it as __main__,
so a plain __main__ guard would never fire -- run_as_freecadcmd_script() at the
bottom is what makes the CLI work, and it is also what stops main() from firing
when something merely IMPORTS this module. That matters now:
src/freecad/build_lead_config.py imports segment_plan() and sweep_segments()
from here so that "sweep a lead along the measured canal centreline" has exactly
ONE implementation in the repo. Before the guard existed, importing this file
built a lead and wrote STLs as a side effect.

freecadcmd also consumes its own argv, hence MAKE_LEAD_ARGS. And it swallows
stdout, so results are written to a report file as well.

For anything beyond a one-off lead -- vertebral levels, lateral offsets, DRG
trajectories, fit validation -- use src/freecad/build_lead_config.py and
src/freecad/lead_defaults.yaml, which drive this module.
"""
import json
import math
import os
import shlex
import sys

import FreeCAD as App
import Part
import Mesh

HERE = os.path.dirname(os.path.abspath(__file__))
MIDLINE_X = 56.60


def polyval(coef, x):
    v = 0.0
    for c in coef:
        v = v * x + c
    return v


def opt(argv, name, cast, default):
    return cast(argv[argv.index(name) + 1]) if name in argv else default


def segment_plan(n_contacts, c_len, gap, tail, z_center):
    """Split a lead into alternating insulator and contact runs along z.

    Returns ([(kind, index, z_start, length), ...], z0, total_length), caudal to
    rostral, where kind is "contact" or "insulator" and index numbers the
    contacts from 1. The two tails are symmetric, so z_center is both the centre
    of the contact array and the centre of the whole lead.

    Segments abut exactly rather than overlapping: neighbours share a flat
    interface, which is what the mesher and a bonded electric contact want, and
    it is how RADO models its own contacts (short full-diameter segments).

    ZERO-LENGTH SEGMENTS ARE DROPPED rather than emitted, because a cylinder of
    zero height is not a shape: Part.makeCylinder() would be handed a length of
    0 and a direction vector it cannot normalise, and the failure comes out as
    an OCC exception several frames away from the `tail: 0.0` that caused it.
    Two leads reach this: one with no tail (RADO's own DRG lead has its
    insulator stop at the end contacts, so reproducing it needs tail 0), and the
    degenerate `gap: 0` where contacts would abut. Every lead that
    existed before this guard has tail 6.0 and gap 1.0, so none of them changes.

    z is the sweep PARAMETER, not necessarily a z coordinate: a DRG lead is
    planned in u, the lateral distance from the midline. The names are kept
    because the rostro-caudal case is the one anyone reads this for.
    """
    total = n_contacts * c_len + (n_contacts - 1) * gap + 2 * tail
    z0 = z_center - 0.5 * total

    segments = []
    z = z0
    if tail > 0:
        segments.append(("insulator", 0, z, tail))
    z += tail
    for i in range(n_contacts):
        segments.append(("contact", i + 1, z, c_len)); z += c_len
        if i < n_contacts - 1:
            if gap > 0:
                segments.append(("insulator", i + 1, z, gap))
            z += gap
    if tail > 0:
        segments.append(("insulator", n_contacts, z, tail))
    return segments, z0, total


def sweep_path(segments, point_of, radius):
    """THE one implementation of "sweep a lead along a measured centreline".

    Each segment becomes a cylinder placed along the LOCAL TANGENT of the
    centreline: from point_of(t0) to point_of(t0 + len). Neighbours therefore
    abut on a flat interface and never overlap, which is what the mesher and a
    bonded electric contact want.

    point_of maps the segment plan's parameter to a point in space, and it is a
    callable because the three lead types run along different axes:

        dorsal / ventral   t is z, and point_of(z) = (x, y(z), z) -- the lead
                           runs rostro-caudally along the epidural canal, whose
                           centre migrates ~7 mm posteriorly over the length of
                           an 8-contact lead because the spine is kyphotic. A
                           lead built at fixed y would walk out of the space and
                           through the dura.
        DRG                t is u, the lateral distance from the midline, and
                           point_of(u) = (midline +/- u, y(u), z(u)) from
                           measure_foramen.py -- the lead runs OUT THROUGH THE
                           FORAMEN, so both of the other two coordinates follow
                           a curve and neither is the sweep parameter.

    Returns [(kind, index, Part.Shape), ...] parallel to `segments`. It builds
    shapes and nothing else -- no document, no files -- so the caller decides
    whether they become a preview in a live document or STLs on disk.
    """
    shapes = []
    for kind, idx, ts, ln in segments:
        a = App.Vector(*point_of(ts))
        b = App.Vector(*point_of(ts + ln))
        d = b.sub(a)
        shapes.append((kind, idx, Part.makeCylinder(radius, d.Length, a, d.normalize())))
    return shapes


def sweep_segments(segments, x, y_of_z, radius):
    """Sweep a rostro-caudal lead: the dorsal/ventral case of sweep_path().

    Kept as its own name because it is the signature every caller of this module
    used before DRG leads existed, and because "at a fixed x, following y(z)" is
    the whole description of an epidural lead. The arithmetic is unchanged --
    the same three vectors in the same order -- so leads generated through here
    are bit-for-bit what they were.

    y_of_z is a callable, deliberately: main() passes the midline polynomial
    from epidural_corridor.json, while build_lead_config.py passes a centreline
    re-measured at the requested lateral offset, where the midline polynomial is
    the wrong curve. Both get identical geometry code.
    """
    return sweep_path(segments, lambda z: (x, y_of_z(z), z), radius)


def fuse_insulator(shapes):
    """Fuse the insulator segments of sweep_segments() output into one shape.

    The contacts stay separate -- each is its own electrode and its own boundary
    condition in the FEM solve -- but the insulator is one physical body, and
    exporting it as six or ten disconnected cylinders would give the mesher
    coincident faces to argue about.

    Raises when a lead has no insulator at all, which segment_plan() will produce
    for `tail: 0` together with `gap: 0`: that is a lead made of bare metal end
    to end, and it should be said so rather than turned into an empty STL.
    """
    parts = [s for kind, _idx, s in shapes if kind == "insulator"]
    if not parts:
        raise ValueError("this lead has no insulator: with tail 0 and gap 0 the "
                         "contacts abut and there is nothing between or beyond "
                         "them. Give it a tail, or a gap between contacts.")
    fused = parts[0]
    for part in parts[1:]:
        fused = fused.fuse(part)
    return fused


def main():
    argv = shlex.split(os.environ.get("MAKE_LEAD_ARGS", ""))
    side = opt(argv, "--side", str, "dorsal")
    n_contacts = opt(argv, "--contacts", int, 8)
    c_len = opt(argv, "--contact-length", float, 3.0)
    gap = opt(argv, "--gap", float, 1.0)
    dia = opt(argv, "--diameter", float, 1.3)
    tail = opt(argv, "--tail", float, 6.0)
    x = opt(argv, "--x", float, MIDLINE_X)
    z_center = opt(argv, "--z-center", float, None)
    tag = opt(argv, "--tag", str, None)

    corridor_path = os.path.join(HERE, "epidural_corridor.json")
    corridor = json.load(open(corridor_path))
    if side not in corridor["sides"]:
        raise SystemExit("side must be dorsal or ventral")
    s = corridor["sides"][side]
    coef = s["centre_poly_coef_high_to_low"]
    z_usable = s["usable_z"]
    if z_center is None:
        z_center = 0.5 * (z_usable[0] + z_usable[1])

    segments, z0, total = segment_plan(n_contacts, c_len, gap, tail, z_center)
    r = 0.5 * dia

    tag = tag or "%s_z%.0f_%dc" % (side, z_center, n_contacts)
    outdir = os.path.join(HERE, "generated_leads", tag)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    lines = []
    log = lines.append
    log("SCS lead: %s, %d contacts x %.2f mm, %.2f mm gaps, %.2f mm dia"
        % (side, n_contacts, c_len, gap, dia))
    log("swept along the measured %s centreline; x = %.2f, z %.2f..%.2f (centre %.2f)"
        % (side, x, z0, z0 + total, z_center))
    log("channel thickness there: min %.2f mm, median %.2f mm"
        % (s["thickness_min_mm"], s["thickness_median_mm"]))
    if dia > s["thickness_min_mm"]:
        log("WARNING: %.2f mm lead exceeds the %.2f mm minimum channel thickness -- it "
            "will overlap dura or the canal wall somewhere along its length."
            % (dia, s["thickness_min_mm"]))
    if z0 < z_usable[0] or z0 + total > z_usable[1]:
        log("WARNING: lead spans z %.1f..%.1f but the centreline is only measured "
            "over %.1f..%.1f; ends are extrapolated."
            % (z0, z0 + total, z_usable[0], z_usable[1]))

    shapes = sweep_segments(segments, x, lambda z: polyval(coef, z), r)

    doc = App.newDocument("scs_lead")
    written = []
    for (kind, idx, shape), (_k, _i, zs, ln) in zip(shapes, segments):
        if kind != "contact":
            continue
        obj = doc.addObject("Part::Feature", "contact%d" % idx)
        obj.Shape = shape
        p = os.path.join(outdir, "SCS Lead Electrode %d.stl" % idx)
        Mesh.export([obj], p)
        written.append((os.path.basename(p), zs, ln, polyval(coef, zs)))

    ins = doc.addObject("Part::Feature", "insulator")
    ins.Shape = fuse_insulator(shapes)
    pi = os.path.join(outdir, "SCS Lead Insulator.stl")
    Mesh.export([ins], pi)
    written.append((os.path.basename(pi), z0, total, polyval(coef, z0)))

    tilt = math.degrees(math.atan2(polyval(coef, z0 + total) - polyval(coef, z0), total))
    log("")
    log("centreline y: %.2f at the caudal end -> %.2f at the rostral end (%.2f mm rise, %.1f deg)"
        % (polyval(coef, z0), polyval(coef, z0 + total),
           polyval(coef, z0 + total) - polyval(coef, z0), tilt))
    log("")
    log("wrote %d files to %s" % (len(written), outdir))
    for name, zs, ln, y in written:
        log("   %-30s z0 %7.2f  len %6.2f  y %6.2f" % (name, zs, ln, y))
    log("")
    log("contact centres (z, y): " + ", ".join(
        "(%.1f, %.2f)" % (sz + 0.5 * sl, polyval(coef, sz + 0.5 * sl))
        for k, _, sz, sl in segments if k == "contact"))

    open(os.path.join(outdir, "lead_report.txt"), "w").write("\n".join(lines) + "\n")
    sys.stdout.write("\n".join(lines) + "\n")


def run_as_freecadcmd_script():
    """True when a FreeCAD interpreter was handed THIS file to run.

    Same guard as apply_colors.py and measure_corridor.py. It exists here so
    that `from make_scs_lead import sweep_segments` does not build a lead and
    write STLs as an import side effect -- see the module docstring.
    """
    if len(sys.argv) < 2:
        return False
    if not os.path.basename(sys.argv[0]).lower().startswith("freecad"):
        return False
    return os.path.abspath(sys.argv[1]) == os.path.abspath(__file__)


if __name__ == "__main__":
    main()
elif run_as_freecadcmd_script():
    main()
