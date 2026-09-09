#!/usr/bin/env python3
"""Colour every body in an NBF_RADO-SCS document by its tissue.

Colours come from `color_rgb` in src/ansys/tissue_map.yaml, which mirrors the
Appearance colour Mohamed set on each material in the Ansys Engineering Data.
The point is that a body looks the same in the FreeCAD tree as it does in
Mechanical, so a mis-assigned body is visible at a glance instead of only
showing up as a wrong conductivity in a solved field.

Which body is which tissue is decided by importing check_tissue_map.py's
parser and its matches() -- not by a second copy of the same logic. If the two
scripts ever disagreed, `check_tissue_map.py` would sign off a mapping that the
picture then contradicts, which is exactly the failure this is meant to catch.
Both are parsed without PyYAML so they run on a bare interpreter; the repo is
public and a collaborator should not have to pip-install anything first.

Transparency does NOT live in the tissue map. The map is a mirror of Engineering
Data -- material names, resistivities, colours -- and Engineering Data has no
notion of transparency, so putting it there would make the map stop being a
faithful mirror. It lives in TRANSPARENCY below instead.


WHICH FREECAD THIS NEEDS
------------------------
Colour and transparency are ViewObject properties, and ViewObjects only exist
when the Gui layer is up. On this build -- FreeCAD 26.3.0 (git 48502), built
from source -- plain `freecadcmd` CANNOT do it:

    import FreeCADGui; FreeCADGui.setupWithoutGUI()

imports and returns cleanly, but leaves FreeCAD.GuiUp == 0, so no view provider
is ever attached and every obj.ViewObject is still None. Verified, not assumed.
Running the full GUI binary offscreen (QT_QPA_PLATFORM=offscreen freecad
script.py) did not come back within two minutes either.

So the write path runs inside a *running* FreeCAD GUI. --dry-run needs no
FreeCAD at all and works under plain python3.

    # plan only, no FreeCAD, writes nothing
    python3 src/freecad/apply_colors.py --dry-run

    # apply, from the Python console of a running FreeCAD (or over the MCP
    # server's execute_code, which is the same interpreter):
    import importlib.util, os
    p = "/home/mohamed/Projects/SCS-Modeling/src/freecad/apply_colors.py"
    os.environ["APPLY_COLORS_ARGS"] = "--doc /home/mohamed/Projects/SCS-Modeling/NBF_RADO-SCS.FCStd --backup"
    spec = importlib.util.spec_from_file_location("apply_colors", p)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    m.main()

    # under freecadcmd it will refuse rather than crash on ViewObject is None:
    freecadcmd src/freecad/apply_colors.py

Options come from the APPLY_COLORS_ARGS environment variable whenever the
interpreter is a FreeCAD one, because FreeCAD owns the command line: it consumes
every flag itself, and passing flags after `--pass` makes freecadcmd skip the
script entirely. Same convention as apply_labels.py.

The other freecadcmd traps apply_labels.py documents also apply here:
freecadcmd IMPORTS the script rather than exec'ing it, so __name__ is never
"__main__" (see run_as_freecadcmd_script at the bottom); print() is swallowed,
so everything goes to a report file (--report) that ends in a "RESULT:" line;
and it segfaults on exit AFTER the document is safely saved, leaving a `core`
file and a meaningless exit code -- read the report, not $?.


THE FREECAD 1.0+ APPEARANCE API
-------------------------------
On this build a Mesh::Feature ViewObject's PropertiesList contains
ShapeAppearance (a tuple of Material) and Transparency, but NOT ShapeColor.
ShapeColor still works as an attribute -- it is a compatibility shim onto
ShapeAppearance[0].DiffuseColor -- and on 26.3.0 a write through it does stick.
apply_appearance() sets ShapeColor, reads the value back off ShapeAppearance,
and only falls back to assigning a Material if the shim did not take, so this
keeps working if the shim is dropped later. The report names which path was
used. Colours are 0.0-1.0 floats in FreeCAD, hence the /255.

Idempotent: it compares the colour and transparency already on each ViewObject
and skips the ones that are right, so a second run reports 0 changed and does
not rewrite the document.
"""
import argparse
import os
import re
import shlex
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))   # src/freecad/ -> repo root
DEFAULT_DOC = os.path.join(REPO, "NBF_RADO-SCS.FCStd")
DEFAULT_STL = os.path.join(REPO, "STL_files")
DEFAULT_MAP = os.path.join(HERE, "..", "ansys", "tissue_map.yaml")
DEFAULT_REPORT = os.path.join(HERE, "apply_colors_report.txt")

# HERE is on sys.path when this file is the script being run, but not when a
# running FreeCAD GUI loads it by path through importlib, which is how the write
# path is invoked -- so put both directories on explicitly.
sys.path.insert(0, os.path.join(HERE, "..", "ansys"))
sys.path.insert(0, HERE)
from check_tissue_map import parse_tissue_map, matches      # noqa: E402
from apply_labels import sanitize                           # noqa: E402

# Per-tissue transparency, 0 = opaque, 100 = invisible. Chosen so the cord reads
# through the shells that wrap it: the eye has to see grey/white matter through
# vertebra + epidural fat + dura + CSF at once, and transparency compounds, so
# the outer shells have to be very transparent for the innermost one to show at
# all. Hardware is opaque because where the contacts sit relative to the cord is
# the whole question. Not in tissue_map.yaml on purpose -- see module docstring.
TRANSPARENCY = {
    "vertebra":             78,
    "intervertebral_disc":  78,
    "soft_tissue":          88,
    "epidural_space":       78,   # the epidural fat the lead sits in
    "meninges_dura":        72,
    "csf":                  72,
    "blood_vessel":         55,   # the thoracic aorta is a big solid cylinder;
                                  # below ~50 it hides everything behind it
    "nerve_root":           58,   # 131 bodies, so they screen the cord in bulk
    "dorsal_root_ganglion": 35,
    "sympathetic_chain":    50,
    "white_matter":         15,
    "grey_matter":           5,
    "electrode_contact":     0,
    "lead_insulation":       0,
}

# Object Names FreeCAD derived from an imported file are turned back into
# something the tissue patterns can be matched against by undoing sanitize():
# every run of non-alphanumerics becomes one space. Only ever used for objects
# the STL index does not already account for -- see lookup().
UNSANITIZE = re.compile(r"[^A-Za-z0-9]+")


def is_freecad_interpreter():
    """True when we are running under FreeCAD rather than a plain python3.

    Covers freecadcmd and the GUI's console, where argv[0] is ".../FreeCAD".
    Lowercased because the GUI binary is "FreeCAD" and the console one is
    "freecadcmd". Two things hang off this: FreeCAD swallows both the command
    line and stdout, so under it options come from an environment variable and
    output goes to the report file instead.
    """
    return os.path.basename(sys.argv[0]).lower().startswith("freecad")


def script_args():
    """Arguments meant for this script, under either interpreter."""
    if is_freecad_interpreter():
        return shlex.split(os.environ.get("APPLY_COLORS_ARGS", ""))
    return sys.argv[1:]


def build_index(stl_dir, tissues):
    """Return (object Name -> (stem, tissue)), plus the ambiguous and unmatched stems.

    Matching is check_tissue_map.py's: every pattern is tried against the
    filename stem, first tissue in file order wins. sanitize() is apply_labels.py's
    reproduction of the object Name FreeCAD derives from that filename.
    """
    stems = sorted(os.path.splitext(f)[0] for f in os.listdir(stl_dir)
                   if f.lower().endswith(".stl"))
    index, ambiguous, unmatched = {}, [], []
    for stem in stems:
        hits = []
        for t in tissues:
            for pat in t["patterns"]:
                if matches(stem, pat):
                    hits.append(t)
                    break
        if not hits:
            unmatched.append(stem)
            continue
        index[sanitize(stem)] = (stem, hits[0])
        if len(hits) > 1:
            ambiguous.append((stem, [t["name"] for t in hits]))
    return index, ambiguous, unmatched, stems


def lookup(name, index, tissues):
    """Resolve one object Name to (stem, tissue, note), or (None, None, note).

    The STL index is tried first, so for every body that came from STL_files the
    answer is byte-for-byte check_tissue_map.py's and the two cannot disagree.

    The fallback matches the tissue patterns against the object Name itself, and
    exists for the 8-contact lead in the dorsal/ventral variants. Those STLs are
    named exactly like RADO's own 4-contact lead files already in the document,
    so FreeCAD uniquifies the Names -- and not the way you would guess: it
    strips the trailing digits off the base name before appending its counter,
    so "SCS Lead Electrode 1.stl" lands as SCS_Lead_Electrode_001, not
    SCS_Lead_Electrode_1001. The contact number in the Name is therefore not
    trustworthy (the Label carries it). What survives intact is the "SCS Lead
    Electrode" part, which is all the pattern needs, so matching the pattern
    against the Name is both simpler and more robust than trying to reverse
    FreeCAD's uniquifier. It uses the same matches() and the same first-in-file
    -wins order as the index path.
    """
    if name in index:
        stem, tissue = index[name]
        return stem, tissue, ""
    key = UNSANITIZE.sub(" ", name).strip()
    for t in tissues:
        for pat in t["patterns"]:
            if matches(key, pat):
                return key, t, "matched on object Name"
    return None, None, "no tissue pattern matched"


def rgb01(color_rgb):
    """0-255 map colour -> the 0.0-1.0 floats FreeCAD wants."""
    return tuple(c / 255.0 for c in color_rgb)


def read_colour(vo):
    """Current diffuse colour as an (r, g, b) triple, read off ShapeAppearance.

    ShapeAppearance is the property that actually exists on this build, so it is
    the honest place to read back from whichever way the value was written.
    """
    return tuple(vo.ShapeAppearance[0].DiffuseColor)[:3]


def same_colour(a, b):
    """Colours equal to within one 0-255 step, which is all the map can express."""
    return all(abs(x - y) < 1.0 / 255.0 for x, y in zip(a, b))


def apply_appearance(vo, rgb, transparency):
    """Set diffuse colour and transparency; return which API the colour went through."""
    used = "ShapeColor"
    try:
        vo.ShapeColor = rgb
        if not same_colour(read_colour(vo), rgb):
            used = None
    except Exception:
        used = None
    if used is None:
        # ShapeColor is only a shim onto ShapeAppearance[0].DiffuseColor; if a
        # future build drops it, write the Material directly. The tuple has to
        # be reassigned -- ShapeAppearance is stored by value, so mutating the
        # Material returned by the getter alone does not stick.
        mat = vo.ShapeAppearance[0]
        mat.DiffuseColor = rgb
        vo.ShapeAppearance = (mat,)
        used = "ShapeAppearance"
    vo.Transparency = int(transparency)
    return used


def open_document(path):
    """Open the document, or reuse it if this FreeCAD already has it open.

    Reopening a file the running GUI already holds would give two in-memory
    copies of a 25 MB document and let one silently overwrite the other's save.
    """
    import FreeCAD
    target = os.path.realpath(path)
    for doc in FreeCAD.listDocuments().values():
        if doc.FileName and os.path.realpath(doc.FileName) == target:
            return doc, True
    return FreeCAD.openDocument(path), False


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--doc", default=DEFAULT_DOC)
    ap.add_argument("--stl-dir", default=DEFAULT_STL)
    ap.add_argument("--map", default=DEFAULT_MAP)
    ap.add_argument("--report", default=DEFAULT_REPORT)
    ap.add_argument("--dry-run", action="store_true",
                    help="report the full plan; touches neither the document nor FreeCAD")
    ap.add_argument("--backup", action="store_true",
                    help="copy the .FCStd to .precolor.bak before writing")
    ap.add_argument("--keep-open", action="store_true",
                    help="leave the document open after saving (for a GUI session)")
    ap.add_argument("--clean-core", action="store_true",
                    help="delete the `core` file the previous freecadcmd run's exit segfault left")
    args = ap.parse_args(script_args())

    out = []
    def say(fmt, *a):
        out.append(fmt % a if a else fmt)

    if args.clean_core and os.path.isfile("core"):
        os.remove("core")
        say("removed the core file left by the previous freecadcmd exit segfault")

    tissues = parse_tissue_map(args.map)
    index, ambiguous, unmatched, stems = build_index(args.stl_dir, tissues)

    say("mode:      %s", "DRY RUN" if args.dry_run else "WRITE")
    say("tissues:   %d from %s", len(tissues), os.path.relpath(args.map, REPO))
    say("STL files: %d in %s", len(stems), os.path.relpath(args.stl_dir, REPO))
    say("document:  %s", args.doc)
    say("")

    # A tissue with a pattern but no colour would silently keep FreeCAD's default
    # grey, which is indistinguishable from a body that was simply missed.
    no_colour = [t["name"] for t in tissues if t["patterns"] and not t["color_rgb"]]
    if no_colour:
        say("ERROR: tissues with patterns but no color_rgb: %s", ", ".join(no_colour))
        say("")
        say("RESULT: aborted, tissue map incomplete")
        write_report(args.report, out)
        return 2

    say("%-24s %-22s %8s %6s  %s", "TISSUE", "ENG-DATA MATERIAL",
        "RGB", "TRANSP", "STL BODIES")
    plan_counts = {}
    for _n, (_stem, t) in index.items():
        plan_counts[t["name"]] = plan_counts.get(t["name"], 0) + 1
    for t in tissues:
        n = plan_counts.get(t["name"], 0)
        rgb = t["color_rgb"] or [0, 0, 0]
        flag = "" if n else "   <-- no bodies matched"
        say("%-24s %-22s %3d,%3d,%3d %5d  %5d%s", t["name"], t["ed_material"] or "-",
            rgb[0], rgb[1], rgb[2], TRANSPARENCY.get(t["name"], 0), n, flag)
    say("%-24s %-22s %8s %6s  %5d", "TOTAL", "", "", "", sum(plan_counts.values()))
    say("")

    missing_transp = [t["name"] for t in tissues if t["name"] not in TRANSPARENCY]
    if missing_transp:
        say("NO TRANSPARENCY DEFAULT (%d) -- would be left opaque: %s",
            len(missing_transp), ", ".join(missing_transp))
        say("")
    if ambiguous:
        say("AMBIGUOUS (%d) -- matched more than one tissue; first in map order wins:",
            len(ambiguous))
        for stem, hits in ambiguous:
            say("   %-58s -> %s", stem[:58], hits)
        say("")
    if unmatched:
        # unmatched_policy: error. A body with no tissue has no colour and no
        # material, so it is a hole in the model, not a cosmetic problem.
        say("UNMATCHED (%d) -- no tissue pattern matched:", len(unmatched))
        for stem in unmatched:
            say("   %s", stem)
        say("")
        say("RESULT: aborted, %d STL bodies unmatched (unmatched_policy: error)",
            len(unmatched))
        write_report(args.report, out)
        return 1

    if args.dry_run:
        say("WOULD SET (stem -> tissue, rgb, transparency), first 25 of %d:", len(index))
        for name in sorted(index)[:25]:
            stem, t = index[name]
            rgb = t["color_rgb"]
            say("   %-58s -> %-22s %3d,%3d,%3d  t=%d", stem[:58], t["name"],
                rgb[0], rgb[1], rgb[2], TRANSPARENCY.get(t["name"], 0))
        say("")
        say("RESULT: dry run OK, %d bodies would be coloured across %d tissues",
            len(index), len(plan_counts))
        write_report(args.report, out)
        return 0

    # -- write path: FreeCAD only imported here so --dry-run works anywhere ---
    import shutil
    import FreeCAD

    if args.backup:
        bak = args.doc + ".precolor.bak"
        n = 1
        while os.path.exists(bak):
            n += 1
            bak = "%s.precolor%d.bak" % (args.doc, n)
        shutil.copy2(args.doc, bak)
        say("backup:    %s", os.path.basename(bak))

    doc, was_open = open_document(args.doc)
    say("opened:    %s, %d objects%s", doc.Name, len(doc.Objects),
        " (already open in this FreeCAD)" if was_open else "")

    # The whole reason this script cannot run under plain freecadcmd. Fail here,
    # with the fix, rather than 245 AttributeErrors on NoneType.
    if any(o.ViewObject is None for o in doc.Objects):
        say("")
        say("ERROR: ViewObject is None -- this interpreter has no Gui layer.")
        say("       FreeCAD.GuiUp = %s. On FreeCAD 26.3.0, FreeCADGui.setupWithoutGUI()",
            FreeCAD.GuiUp)
        say("       returns cleanly but does not attach view providers, so colours")
        say("       cannot be set from freecadcmd. Run this from the Python console")
        say("       of a running FreeCAD GUI -- see the module docstring.")
        say("")
        say("RESULT: aborted, no Gui layer")
        write_report(args.report, out)
        if not was_open:
            FreeCAD.closeDocument(doc.Name)
        return 3

    changed, already, apis, no_tissue, samples = 0, 0, set(), [], []
    applied_counts, name_matched = {}, []
    for obj in doc.Objects:
        vo = obj.ViewObject
        if vo is None or "ShapeAppearance" not in vo.PropertiesList:
            continue                      # not a drawable body (none in this doc)
        stem, tissue, note = lookup(obj.Name, index, tissues)
        if tissue is None:
            no_tissue.append((obj.Name, obj.Label))
            continue
        if note:
            name_matched.append((obj.Name, tissue["name"]))
        applied_counts[tissue["name"]] = applied_counts.get(tissue["name"], 0) + 1

        rgb = rgb01(tissue["color_rgb"])
        transp = TRANSPARENCY.get(tissue["name"], 0)
        if same_colour(read_colour(vo), rgb) and vo.Transparency == transp:
            already += 1
            continue
        apis.add(apply_appearance(vo, rgb, transp))
        changed += 1
        if len(samples) < 20:
            samples.append((obj.Label, tissue["name"], tissue["color_rgb"], transp))

    say("")
    say("%-24s %s", "TISSUE", "BODIES COLOURED")
    for t in tissues:
        n = applied_counts.get(t["name"], 0)
        if n:
            say("%-24s %5d", t["name"], n)
    say("%-24s %5d", "TOTAL", sum(applied_counts.values()))
    say("")
    say("changed:   %d", changed)
    say("unchanged: %d (already the right colour and transparency)", already)
    say("api used:  %s", ", ".join(sorted(apis)) if apis else "-")
    if name_matched:
        say("")
        say("RESOLVED ON OBJECT NAME (%d) -- no STL of that name in --stl-dir,",
            len(name_matched))
        say("so the tissue patterns were matched against the Name instead:")
        for name, tname in name_matched:
            say("   %-40s -> %s", name, tname)
    if no_tissue:
        say("")
        say("OBJECTS WITH NO TISSUE (%d) -- left alone:", len(no_tissue))
        for name, label in no_tissue:
            say("   %-40s Label=%r", name[:40], label)

    if samples:
        say("")
        say("CHANGED (first %d):", len(samples))
        for label, tname, rgb, transp in samples:
            say("   %-46s %-22s %3d,%3d,%3d  t=%d", label[:46], tname,
                rgb[0], rgb[1], rgb[2], transp)

    if changed:
        doc.save()
        say("")
        say("saved:     %s", args.doc)
    else:
        say("")
        say("no changes, document not rewritten")
    if not args.keep_open and not was_open:
        FreeCAD.closeDocument(doc.Name)

    say("")
    say("RESULT: %d coloured, %d already correct, %d objects without a tissue",
        changed, already, len(no_tissue))
    write_report(args.report, out)
    return 0 if not no_tissue else 1


def write_report(path, lines):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    # Under freecadcmd and inside the GUI, print() goes nowhere, which is why the
    # report file exists at all. Under a plain interpreter it does work, so echo
    # there rather than leaving `python3 apply_colors.py --dry-run` looking like
    # it silently did nothing.
    if not is_freecad_interpreter():
        sys.stdout.write("\n".join(lines) + "\n")


def run_as_freecadcmd_script():
    """True when freecadcmd was handed THIS file to run.

    freecadcmd imports a script rather than exec'ing it, so __name__ is
    "apply_colors" and the usual __main__ guard never fires -- the script would
    load, define everything, and silently do nothing. Comparing the path
    freecadcmd was given against __file__ is what makes
    `freecadcmd src/freecad/apply_colors.py` actually run, without also firing
    when something else imports this module.
    """
    if len(sys.argv) < 2 or not is_freecad_interpreter():
        return False
    return os.path.abspath(sys.argv[1]) == os.path.abspath(__file__)


if __name__ == "__main__":
    sys.exit(main())
elif run_as_freecadcmd_script():
    # No sys.exit here: raising SystemExit out of an imported module makes
    # freecadcmd print a traceback over a run that actually succeeded.
    main()
