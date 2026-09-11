#!/usr/bin/env python3
"""Put every body of an NBF_RADO-SCS document into a group named for its tissue.

This is FreeCAD's answer to Ansys' named selections. Mohamed's actual complaint
was that hiding "all the vertebrae" or "all the vasculature" to get a clear look
at the lead and the cord meant clicking spacebar 123 times. An
App::DocumentObjectGroup is a real container in the tree: select it, press
spacebar, and every body inside it hides at once.

Groups are a good fit here specifically because tissue membership is a clean
PARTITION. An App::DocumentObjectGroup has exclusive parentage -- an object can
live in exactly one -- which would make groups useless for overlapping
selections (say "everything near the electrodes"). But src/ansys/tissue_map.yaml
maps every body to exactly one tissue, and src/ansys/check_tissue_map.py fails
the build if any body matches zero or two entries, so there is nothing to
overlap. Anything that needs overlapping sets wants the macro instead --
src/freecad/tissue_visibility.FCMacro, which resolves bodies by tissue on the
fly and does not care about the tree at all.

Which body is which tissue is decided by importing apply_colors.py's
build_index() and lookup(), which in turn import check_tissue_map.py's parser
and matches(). One resolver, three consumers: the Ansys material check, the
colours, and now the groups. If a body's colour and its group ever disagreed,
the tree would be lying about the physics, and that is the whole failure this
chain of imports exists to make impossible.

WHY THE LABELS ARE NUMBERED
---------------------------
Group Labels are "01 Vertebrae", "02 Discs", ... "11 Vasculature" rather
than plain names. FreeCAD's tree can be switched to alphabetical sort (right-
click the document -> Tree view options), and several people leave it that way;
without a prefix that sort scatters the anatomy -- CSF, Discs, DRG, Dura, Grey
matter -- so you cannot read the model outward from bone to cord. The two-digit
prefix makes the alphabetical order and the anatomical order the same one:
outside in (vertebrae, discs, epidural space, dura, CSF, white, grey), then the
peripheral structures, then the hardware last so the lead is always at the
bottom of the tree where you can find it. Two digits, not one, because there are
more than nine groups and "10" would otherwise sort between "1" and "2".

A NUMBER WITH NO GROUP IS NORMAL
--------------------------------
The numbering comes from position in GROUP_ORDER, not from how many groups there
are, so it has gaps: a tissue with no body in the document gets no group, and
the next tissue keeps its own number rather than shuffling up. In the anatomy
document today that is 12 (soft tissue, whose bodies are not in this model),
13 (lead contacts) and 14 (lead insulation) -- the document holds ANATOMY ONLY
and every lead is a disposable preview, so it ends at "11 Vasculature" with
240 bodies in 11 groups. Existing empty groups are REMOVED on a re-run for the
same reason. The lead tissues stay in tissue_map.yaml regardless; previews and
exported STLs need their colours, and the preview puts its own bodies in a
group of its own ("13 SCS Leads") which is not one of these.


WHICH FREECAD THIS NEEDS -- READ THIS BEFORE REACHING FOR freecadcmd
--------------------------------------------------------------------
Grouping itself is a pure App-level operation: addObject("App::DocumentObject-
Group") and group.addObject(body) both work perfectly under plain freecadcmd,
with FreeCAD.GuiUp == 0 and every ViewObject None. Verified on this build,
FreeCAD 26.3.0 (git 48502).

The SAVE is what you cannot do headless. A document saved by freecadcmd is
written with NO GuiDocument.xml and none of the ShapeAppearance blobs at all --
measured on a copy of the model document:

    before   Document.xml + GuiDocument.xml (894 839 bytes) + ShapeAppearance*0-253
    after    Document.xml only, 255 entries -> 1

There is no Gui layer to serialise, so FreeCAD silently writes the App half of
the file and drops the rest. Every colour and transparency apply_colors.py set
would be gone -- and the document would still open, still contain every body,
and look like a fresh grey import. That is a much worse outcome than a crash, so
this script REFUSES to save when FreeCAD.GuiUp is false rather than producing
it. --force-headless-save exists only for a document that has no view data worth
keeping.

drop_rado_lead.py, the only other script here that saves, refuses on the same
test and for the same reason.

So, like apply_colors.py, the write path runs inside a *running* FreeCAD GUI.
--dry-run needs no FreeCAD at all and works under plain python3.

    # plan only, no FreeCAD, writes nothing
    python3 src/freecad/make_tissue_groups.py --dry-run

    # apply, from the Python console of a running FreeCAD (or over the MCP
    # server's execute_code, which is the same interpreter):
    import importlib.util, os
    p = "/home/mohamed/Projects/SCS-Modeling/src/freecad/make_tissue_groups.py"
    os.environ["MAKE_TISSUE_GROUPS_ARGS"] = (
        "--doc /home/mohamed/Projects/SCS-Modeling/NBF_RADO-SCS.FCStd --backup")
    spec = importlib.util.spec_from_file_location("make_tissue_groups", p)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    m.main()

Options come from the MAKE_TISSUE_GROUPS_ARGS environment variable under any
FreeCAD interpreter, because FreeCAD owns the command line: it consumes every
flag itself, and passing flags after `--pass` makes freecadcmd skip the script
entirely. Same convention as apply_colors.py and apply_labels.py. Output goes to
a report file ending in a "RESULT:" line because FreeCAD swallows print().

Idempotent: a group is reused if it already exists (matched on the object Name,
which is stable, not the Label, which a user may retitle), and a body already
sitting in the right group is left untouched. A second run reports 0 moved and
does not rewrite the document. Regrouping does not touch the bodies themselves:
Label, Placement, Visibility, ShapeAppearance and Transparency are all
properties OF the body, and membership of a group is a property of the GROUP
(its Group link list), so moving one cannot perturb the other. check_untouched()
below asserts that rather than trusting it.
"""
import argparse
import os
import shlex
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))   # src/freecad/ -> repo root
DEFAULT_DOC = os.path.join(REPO, "NBF_RADO-SCS.FCStd")
DEFAULT_STL = os.path.join(REPO, "STL_files")
DEFAULT_MAP = os.path.join(HERE, "..", "ansys", "tissue_map.yaml")
DEFAULT_REPORT = os.path.join(HERE, "make_tissue_groups_report.txt")

# HERE is on sys.path when this file is the script being run, but not when a
# running FreeCAD GUI loads it by path through importlib, which is how the write
# path is invoked -- so put it on explicitly. apply_colors.py does the same for
# ../ansys, and importing it is what puts check_tissue_map on the path too.
sys.path.insert(0, HERE)
from apply_colors import (build_index, lookup,                     # noqa: E402
                          is_freecad_interpreter, open_document)
from check_tissue_map import parse_tissue_map                      # noqa: E402

# tissue name in tissue_map.yaml -> the human-readable half of the group Label.
# The order of this list IS the numbering: outside in, then peripheral, then
# hardware. See "WHY THE LABELS ARE NUMBERED" above. A tissue in the map but not
# here is still grouped -- it lands at the end with the next number -- so adding
# a tissue to the map cannot silently drop its bodies out of the tree.
GROUP_ORDER = [
    ("vertebra",             "Vertebrae"),
    ("intervertebral_disc",  "Discs"),
    ("epidural_space",       "Epidural space"),
    ("meninges_dura",        "Dura and meninges"),
    ("csf",                  "CSF"),
    ("white_matter",         "White matter"),
    ("grey_matter",          "Grey matter"),
    ("nerve_root",           "Nerve roots"),
    ("dorsal_root_ganglion", "DRG"),
    ("sympathetic_chain",    "Sympathetic chain"),
    ("blood_vessel",         "Vasculature"),
    ("soft_tissue",          "Soft tissue"),
    ("electrode_contact",    "Lead contacts"),
    ("lead_insulation",      "Lead insulation"),
]

# Prefix for the group objects' internal Names. The Name, not the Label, is what
# makes a re-run idempotent: FreeCAD guarantees Names are unique and immutable
# for the life of the object, whereas a Label is free text Mohamed may well
# retitle. Deliberately not a tissue Label prefix -- it never shows in the tree.
NAME_PREFIX = "TissueGroup_"

# The body properties regrouping must not perturb. Read before and after the
# moves and compared exactly; see check_untouched(). Placement is compared via
# its string form because two Placement objects do not compare equal by value.
WATCHED = ("Label", "Visibility", "Placement")


def script_args():
    """Arguments meant for this script; see apply_colors.script_args."""
    if is_freecad_interpreter():
        return shlex.split(os.environ.get("MAKE_TISSUE_GROUPS_ARGS", ""))
    return sys.argv[1:]


def group_name_for(tissue):
    """Internal object Name for one tissue's group -- stable, unique, immutable."""
    return NAME_PREFIX + tissue


def group_plan(tissues):
    """Return [(tissue name, group Name, group Label), ...] in tree order.

    Numbering is positional, so it is derived here rather than written out by
    hand: a tissue added to GROUP_ORDER renumbers everything after it in one
    place. Tissues present in the map but missing from GROUP_ORDER are appended
    in map order so they still get a group.
    """
    ordered = [name for name, _ in GROUP_ORDER]
    extra = [t["name"] for t in tissues if t["name"] not in ordered]
    titles = dict(GROUP_ORDER)
    plan = []
    for i, tissue in enumerate(ordered + extra, start=1):
        title = titles.get(tissue) or tissue.replace("_", " ").capitalize()
        plan.append((tissue, group_name_for(tissue), "%02d %s" % (i, title)))
    return plan


def snapshot(objs):
    """{Name: (Label, Visibility, str(Placement))} for the untouched-check."""
    shot = {}
    for obj in objs:
        shot[obj.Name] = tuple(
            str(getattr(obj, p)) if p == "Placement" else getattr(obj, p)
            for p in WATCHED)
    return shot


def check_untouched(before, after):
    """Return [(Name, property, before, after), ...] for anything regrouping changed.

    Regrouping SHOULD be a no-op on the bodies -- group membership is a link
    held by the group, not a property of the member -- but "should" is how a
    document quietly loses 245 Labels. Cheap to check, so it is checked.
    """
    diffs = []
    for name, vals in before.items():
        if name not in after:
            diffs.append((name, "-", "present", "GONE"))
            continue
        for prop, was, now in zip(WATCHED, vals, after[name]):
            if was != now:
                diffs.append((name, prop, was, now))
    return diffs


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
                    help="copy the .FCStd to .pregroup.bak before writing")
    ap.add_argument("--keep-open", action="store_true",
                    help="leave the document open after saving (for a GUI session)")
    ap.add_argument("--force-headless-save", action="store_true",
                    help="save even with no Gui layer -- DISCARDS every colour, "
                         "transparency and visibility flag in the document")
    args = ap.parse_args(script_args())

    out = []
    def say(fmt, *a):
        out.append(fmt % a if a else fmt)

    tissues = parse_tissue_map(args.map)
    index, ambiguous, unmatched, stems = build_index(args.stl_dir, tissues)
    plan = group_plan(tissues)

    say("mode:      %s", "DRY RUN" if args.dry_run else "WRITE")
    say("tissues:   %d from %s", len(tissues), os.path.relpath(args.map, REPO))
    say("STL files: %d in %s", len(stems), os.path.relpath(args.stl_dir, REPO))
    say("document:  %s", args.doc)
    say("")

    if unmatched:
        # Same policy as apply_colors.py: a body with no tissue has no group and
        # no material, so it is a hole in the model, not a cosmetic problem.
        say("UNMATCHED (%d) -- no tissue pattern matched:", len(unmatched))
        for stem in unmatched:
            say("   %s", stem)
        say("")
        say("RESULT: aborted, %d STL bodies unmatched (unmatched_policy: error)",
            len(unmatched))
        write_report(args.report, out)
        return 1

    stl_counts = {}
    for _n, (_stem, t) in index.items():
        stl_counts[t["name"]] = stl_counts.get(t["name"], 0) + 1

    say("Counted from the STL FILES, which is all a dry run can see. The anatomy")
    say("document holds fewer: drop_rado_lead.py took RADO's five lead bodies out")
    say("of it, so 13 and 14 below have no bodies THERE and a real run neither")
    say("creates nor keeps them. The STLs stay on disk either way.")
    say("")
    say("%-24s %-32s %s", "GROUP LABEL", "OBJECT NAME", "STL BODIES")
    for tissue, gname, glabel in plan:
        n = stl_counts.get(tissue, 0)
        flag = "" if n else "   <-- no bodies: not created, and removed if present"
        say("%-24s %-32s %5d%s", glabel, gname, n, flag)
    say("%-24s %-32s %5d", "TOTAL", "", sum(stl_counts.values()))
    say("")

    if ambiguous:
        say("AMBIGUOUS (%d) -- matched more than one tissue; first in map order wins,",
            len(ambiguous))
        say("which is also the group each lands in:")
        for stem, hits in ambiguous:
            say("   %-58s -> %s", stem[:58], hits)
        say("")

    if args.dry_run:
        say("RESULT: dry run OK, %d groups over %d bodies",
            len([1 for t, _g, _l in plan if stl_counts.get(t)]), sum(stl_counts.values()))
        write_report(args.report, out)
        return 0

    # -- write path: FreeCAD only imported here so --dry-run works anywhere ---
    import shutil
    import FreeCAD

    # Refuse BEFORE touching anything. A headless save would silently strip
    # GuiDocument.xml and every ShapeAppearance blob out of the document -- see
    # the module docstring. Checked here rather than after the moves so a
    # mistaken freecadcmd invocation costs nothing.
    if not FreeCAD.GuiUp and not args.force_headless_save:
        say("ERROR: FreeCAD.GuiUp = %s -- no Gui layer.", FreeCAD.GuiUp)
        say("       The grouping itself would work fine headless, but the SAVE would")
        say("       write Document.xml alone: no GuiDocument.xml, no ShapeAppearance")
        say("       blobs. Every colour and transparency apply_colors.py set, and the")
        say("       Visibility flags on RADO's 4-contact lead, would be lost -- and the")
        say("       document would still open and look like a fresh grey import.")
        say("       Run this from the Python console of a running FreeCAD GUI instead;")
        say("       see the module docstring. --force-headless-save overrides, and is")
        say("       only correct on a document with no view data worth keeping.")
        say("")
        say("RESULT: aborted, no Gui layer (would destroy the document's view data)")
        write_report(args.report, out)
        return 3

    if args.backup:
        bak = args.doc + ".pregroup.bak"
        n = 1
        while os.path.exists(bak):
            n += 1
            bak = "%s.pregroup%d.bak" % (args.doc, n)
        shutil.copy2(args.doc, bak)
        say("backup:    %s", os.path.basename(bak))

    doc, was_open = open_document(args.doc)
    say("opened:    %s, %d objects%s", doc.Name, len(doc.Objects),
        " (already open in this FreeCAD)" if was_open else "")

    bodies = [o for o in doc.Objects if o.TypeId != "App::DocumentObjectGroup"]
    before = snapshot(bodies)

    # Which group each body should end up in. Resolution is apply_colors.py's,
    # so a body's group and its colour are decided by the same call.
    wanted, no_tissue = {}, []
    for obj in bodies:
        _stem, tissue, _note = lookup(obj.Name, index, tissues)
        if tissue is None:
            no_tissue.append((obj.Name, obj.Label))
            continue
        wanted[obj.Name] = tissue["name"]

    counts = {}
    for tname in wanted.values():
        counts[tname] = counts.get(tname, 0) + 1

    # Create the groups in plan order so the tree's *creation* order matches the
    # numbering too, for anyone who has not switched the tree to alphabetical.
    groups, created, reused, relabelled = {}, [], [], []
    for tissue, gname, glabel in plan:
        if not counts.get(tissue):
            continue
        g = doc.getObject(gname)
        if g is None:
            g = doc.addObject("App::DocumentObjectGroup", gname)
            created.append(glabel)
        else:
            reused.append(glabel)
        if g.Label != glabel:
            if g.Label not in ("", gname):
                relabelled.append((g.Label, glabel))
            g.Label = glabel
        groups[tissue] = g

    # Move. An App::DocumentObjectGroup has exclusive parentage, so a body that
    # is already in some OTHER group has to be removed from it first -- FreeCAD
    # will otherwise happily list it under both and the tree grows a duplicate.
    moved, already = 0, 0
    for obj in bodies:
        tname = wanted.get(obj.Name)
        if tname is None:
            continue
        g = groups[tname]
        parent = obj.getParentGroup() if hasattr(obj, "getParentGroup") else None
        if parent is not None and parent.Name == g.Name:
            already += 1
            continue
        if parent is not None:
            parent.removeObject(obj)
        g.addObject(obj)
        moved += 1

    # An EMPTY group in the plan is removed rather than left standing. The
    # anatomy document holds no lead any more (drop_rado_lead.py), so
    # "13 Lead contacts" and "14 Lead insulation" would otherwise sit in the
    # tree announcing hardware that is not in the document -- the same lie the
    # numbering exists to prevent. Only groups this script owns (NAME_PREFIX)
    # and only when they hold nothing, so a group someone put a body in by hand
    # survives. The tissues themselves stay in tissue_map.yaml: a previewed or
    # exported lead still needs their colours.
    emptied = []
    for tissue, gname, glabel in plan:
        g = doc.getObject(gname)
        if g is not None and not g.Group and not counts.get(tissue):
            emptied.append(glabel)
            doc.removeObject(gname)
            groups.pop(tissue, None)

    doc.recompute()
    after = snapshot(bodies)
    diffs = check_untouched(before, after)

    say("")
    say("%-24s %-32s %s", "GROUP LABEL", "OBJECT NAME", "BODIES")
    total = 0
    for tissue, gname, glabel in plan:
        if tissue not in groups:
            continue
        n = len(groups[tissue].Group)
        total += n
        say("%-24s %-32s %5d", glabel, gname, n)
    say("%-24s %-32s %5d", "TOTAL", "", total)
    say("")
    say("groups created: %d  (%s)", len(created), ", ".join(created) or "-")
    say("groups reused:  %d  (%s)", len(reused), ", ".join(reused) or "-")
    say("groups removed: %d  (%s)   empty, no body of that tissue in the document",
        len(emptied), ", ".join(emptied) or "-")
    if relabelled:
        say("groups renamed: %s",
            ", ".join("%r -> %r" % (a, b) for a, b in relabelled))
    say("bodies moved:   %d", moved)
    say("bodies already in the right group: %d", already)

    if no_tissue:
        say("")
        say("OBJECTS WITH NO TISSUE (%d) -- left ungrouped at the document root:",
            len(no_tissue))
        for name, label in no_tissue:
            say("   %-40s Label=%r", name[:40], label)

    say("")
    if diffs:
        say("REGROUPING CHANGED BODY PROPERTIES (%d) -- this should be empty:", len(diffs))
        for name, prop, was, now in diffs[:40]:
            say("   %-32s %-12s %r -> %r", name[:32], prop, was, now)
    else:
        say("untouched-check: %d bodies, %s all unchanged by the regrouping",
            len(before), "/".join(WATCHED))

    if moved or created or relabelled:
        doc.save()
        say("saved:     %s", args.doc)
    else:
        say("no changes, document not rewritten")
    if not args.keep_open and not was_open:
        FreeCAD.closeDocument(doc.Name)

    say("")
    say("RESULT: %d groups, %d bodies moved, %d already grouped, %d ungrouped, "
        "%d property changes", len(groups), moved, already, len(no_tissue), len(diffs))
    write_report(args.report, out)
    return 0 if not (no_tissue or diffs) else 1


def write_report(path, lines):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    # Under freecadcmd and inside the GUI, print() goes nowhere, which is why the
    # report file exists at all. Under a plain interpreter it does work, so echo
    # there rather than leaving `--dry-run` looking like it silently did nothing.
    if not is_freecad_interpreter():
        sys.stdout.write("\n".join(lines) + "\n")


def run_as_freecadcmd_script():
    """True when a FreeCAD interpreter was handed THIS file to run.

    freecadcmd imports a script rather than exec'ing it, so __name__ is
    "make_tissue_groups" and the usual __main__ guard never fires -- the script
    would load, define everything, and silently do nothing. See apply_colors.py.
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
