#!/usr/bin/env python3
"""Generate a parametric percutaneous SCS lead as STL bodies.

Why this exists
---------------
The lead that ships with RADO is not the one we want to study. Measured from the
geometry (ansys/check_laterality.py), it sits 11-25 mm LEFT of midline at
z = 93-97 mm, draped over the left DRG column -- a DRG lead, not a dorsal-column
lead -- and it has only 4 contacts on a curved insulator. Khadka et al. Figure 1
shows the clinical arrangement instead: a straight percutaneous lead with 8
cylindrical contacts, placed in the posterior epidural space at midline.

This builds that lead, parametrically, so it can be dropped in at any
rostro-caudal level and on either the dorsal or the ventral side -- which is
exactly the sweep the dorsal-vs-ventral study needs.

Anatomical frame (derived from the geometry, NOT from Khadka's filename tags,
which are inverted for most bodies -- see ansys/check_laterality.py):

    +X  anatomical LEFT      midline at x = 56.60 mm
    +Y  posterior / DORSAL   (-Y is anterior/ventral)
    +Z  rostral

Measured corridors at midline (from the RADO STLs):

    dural sac (meninges)      y = 59.0 .. 87.6
    epidural space            y = 57.2 .. 89.9
    => dorsal epidural gap    y = 87.6 .. 89.9   (2.3 mm)
    => ventral epidural gap   y = 57.2 .. 59.0   (1.8 mm)
    usable rostro-caudal span z = 59.4 .. 164.9  (105.5 mm)

The lead axis runs along Z. The body is split along its length into alternating
insulator and contact segments, each a full-diameter cylinder, so the parts
share flat interfaces and never overlap -- which is what the mesher and the
bonded electric contact want. That also matches how the existing RADO contacts
are modelled (short full-diameter segments, ~1.3 mm).

Defaults follow a clinical percutaneous lead (8 contacts, 3 mm long, 1 mm gaps,
1.3 mm diameter). Everything is overridable.

Usage (needs FreeCAD; run headless):

    freecadcmd ansys/make_scs_lead.py            # dorsal, centred, default 8 contacts
    MAKE_LEAD_ARGS="--side ventral --z-center 110" freecadcmd ansys/make_scs_lead.py
    MAKE_LEAD_ARGS="--help" freecadcmd ansys/make_scs_lead.py

NOTE: freecadcmd *imports* this file rather than running it as __main__, so a
`if __name__ == "__main__"` guard would never fire -- the work is done at import
time. It also consumes its own argv, so options come from the environment
variable MAKE_LEAD_ARGS instead of the command line. (Both traps are documented
in docs/ansys_hpc_training_notes.md.) Output is written to a report file because
freecadcmd swallows stdout.
"""
import argparse
import os
import shlex
import sys

import FreeCAD as App
import Part
import Mesh

# Frame constants measured from the RADO STL set (mm).
MIDLINE_X = 56.60
DORSAL_GAP = (87.6, 89.9)     # dura outer .. epidural outer, at midline
VENTRAL_GAP = (57.2, 59.0)
Z_SPAN = (59.4, 164.9)


def parse_args():
    ap = argparse.ArgumentParser(prog="make_scs_lead.py", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--contacts", type=int, default=8, help="number of contacts (default 8)")
    ap.add_argument("--contact-length", type=float, default=3.0, help="mm (default 3.0)")
    ap.add_argument("--gap", type=float, default=1.0, help="mm between contacts (default 1.0)")
    ap.add_argument("--diameter", type=float, default=1.3, help="mm (default 1.3)")
    ap.add_argument("--tail", type=float, default=6.0,
                    help="mm of plain insulator beyond the end contacts (default 6.0)")
    ap.add_argument("--side", choices=("dorsal", "ventral"), default="dorsal")
    ap.add_argument("--x", type=float, default=None, help="mm; default = midline")
    ap.add_argument("--y", type=float, default=None,
                    help="mm; default = centre of the chosen epidural gap")
    ap.add_argument("--z-center", type=float, default=None,
                    help="mm; default = centre of the usable span")
    ap.add_argument("--outdir", default=None, help="default ansys/generated_leads/<tag>")
    ap.add_argument("--tag", default=None, help="name for this configuration")
    return ap.parse_args(shlex.split(os.environ.get("MAKE_LEAD_ARGS", "")))


def main():
    args = parse_args()
    here = os.path.dirname(os.path.abspath(__file__))

    gap_lo, gap_hi = DORSAL_GAP if args.side == "dorsal" else VENTRAL_GAP
    x = MIDLINE_X if args.x is None else args.x
    y = 0.5 * (gap_lo + gap_hi) if args.y is None else args.y
    zc = 0.5 * (Z_SPAN[0] + Z_SPAN[1]) if args.z_center is None else args.z_center

    pitch = args.contact_length + args.gap
    active = args.contacts * args.contact_length + (args.contacts - 1) * args.gap
    total = active + 2 * args.tail
    z0 = zc - 0.5 * total
    r = 0.5 * args.diameter

    tag = args.tag or "%s_z%.0f_%dc" % (args.side, zc, args.contacts)
    outdir = args.outdir or os.path.join(here, "generated_leads", tag)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    lines = []
    def log(m):
        lines.append(str(m))

    log("SCS lead: %s, %d contacts, %.2f mm long, %.2f mm gaps, %.2f mm dia"
        % (args.side, args.contacts, args.contact_length, args.gap, args.diameter))
    log("axis at x=%.2f y=%.2f (gap %.1f..%.1f), z %.2f..%.2f (centre %.2f)"
        % (x, y, gap_lo, gap_hi, z0, z0 + total, zc))
    if not (gap_lo <= y - r and y + r <= gap_hi):
        log("WARNING: a %.2f mm lead does not fit the %.2f mm %s gap -- it will "
            "overlap dura or the epidural boundary and need a boolean, or a "
            "smaller --diameter." % (args.diameter, gap_hi - gap_lo, args.side))
    if z0 < Z_SPAN[0] or z0 + total > Z_SPAN[1]:
        log("WARNING: the lead runs outside the modelled span z=%.1f..%.1f"
            % Z_SPAN)

    # Alternating segments along Z: tail, [contact, gap] x N-1, contact, tail.
    # Each is a full-diameter cylinder so neighbours share a flat interface.
    segments = []            # (kind, index, z_start, length)
    z = z0
    segments.append(("insulator", 0, z, args.tail)); z += args.tail
    for i in range(args.contacts):
        segments.append(("contact", i + 1, z, args.contact_length))
        z += args.contact_length
        if i < args.contacts - 1:
            segments.append(("insulator", i + 1, z, args.gap)); z += args.gap
    segments.append(("insulator", args.contacts, z, args.tail)); z += args.tail

    doc = App.newDocument("scs_lead")
    insulator_parts, written = [], []
    for kind, idx, zs, ln in segments:
        cyl = Part.makeCylinder(r, ln, App.Vector(x, y, zs), App.Vector(0, 0, 1))
        if kind == "contact":
            p = os.path.join(outdir, "SCS Lead Electrode %d.stl" % idx)
            obj = doc.addObject("Part::Feature", "contact%d" % idx)
            obj.Shape = cyl
            Mesh.export([obj], p)
            written.append((p, ln, zs))
        else:
            insulator_parts.append(cyl)

    fused = insulator_parts[0]
    for s in insulator_parts[1:]:
        fused = fused.fuse(s)
    ins = doc.addObject("Part::Feature", "insulator")
    ins.Shape = fused
    pi = os.path.join(outdir, "SCS Lead Insulator.stl")
    Mesh.export([ins], pi)
    written.append((pi, total - args.contacts * args.contact_length, z0))

    log("")
    log("wrote %d files to %s" % (len(written), outdir))
    for p, ln, zs in written:
        log("   %-34s  len %6.2f mm  z0 %7.2f  (%d bytes)"
            % (os.path.basename(p), ln, zs, os.path.getsize(p)))
    log("")
    log("contact centres (z, mm): " +
        ", ".join("%.2f" % (s[2] + 0.5 * s[3]) for s in segments if s[0] == "contact"))

    rep = os.path.join(outdir, "lead_report.txt")
    open(rep, "w").write("\n".join(lines) + "\n")
    sys.stdout.write("\n".join(lines) + "\n")


main()
