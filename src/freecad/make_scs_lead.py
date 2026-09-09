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
so the work happens at import time and a __main__ guard would never fire. It
also consumes its own argv, hence MAKE_LEAD_ARGS. And it swallows stdout, so
results are written to a report file as well.
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

    pitch_total = n_contacts * c_len + (n_contacts - 1) * gap
    total = pitch_total + 2 * tail
    z0 = z_center - 0.5 * total
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

    segments = []
    z = z0
    segments.append(("insulator", 0, z, tail)); z += tail
    for i in range(n_contacts):
        segments.append(("contact", i + 1, z, c_len)); z += c_len
        if i < n_contacts - 1:
            segments.append(("insulator", i + 1, z, gap)); z += gap
    segments.append(("insulator", n_contacts, z, tail))

    doc = App.newDocument("scs_lead")
    insulator_parts, written = [], []
    for kind, idx, zs, ln in segments:
        # place each segment along the local tangent of the centreline
        a = App.Vector(x, polyval(coef, zs), zs)
        b = App.Vector(x, polyval(coef, zs + ln), zs + ln)
        d = b.sub(a)
        cyl = Part.makeCylinder(r, d.Length, a, d.normalize())
        if kind == "contact":
            obj = doc.addObject("Part::Feature", "contact%d" % idx)
            obj.Shape = cyl
            p = os.path.join(outdir, "SCS Lead Electrode %d.stl" % idx)
            Mesh.export([obj], p)
            written.append((os.path.basename(p), zs, ln, a.y))
        else:
            insulator_parts.append(cyl)

    fused = insulator_parts[0]
    for part in insulator_parts[1:]:
        fused = fused.fuse(part)
    ins = doc.addObject("Part::Feature", "insulator")
    ins.Shape = fused
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


main()
