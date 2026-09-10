#!/usr/bin/env python3
"""Remove RADO's own 4-contact DRG lead from NBF_RADO-SCS.FCStd. GUI only.

WHY THE ANATOMY DOCUMENT SHOULD NOT CONTAIN A LEAD
--------------------------------------------------
There is now one model document and it holds ANATOMY ONLY. Every lead is
configuration (lead_configs.yaml) plus geometry generated on demand, previewed
into the open document and thrown away. RADO's five hardware bodies -- four
platinum contacts and an insulator sheath -- are the last thing in the document
that does not belong to that scheme, and leaving them there is not merely untidy:

    Simpleware or Gmsh will mesh them. The FEM solve will treat four metal
    cylinders and a sheath as conductors sitting in the field a few millimetres
    off the cord, and the potentials -- and every axon threshold computed from
    them -- will be wrong.

That is the same reasoning that deleted them from the old dorsal/ventral study
variants, applied to the document those variants were copies of. Hiding a body
(Visibility = False) changes none of it; it only stops FreeCAD drawing it.

WHAT IS AND IS NOT LOST
-----------------------
Nothing is lost. Those five bodies are imports of five STLs that RADO ships and
that this repo tracks unmodified:

    STL_files/SCS Lead Electrode 1.stl  ..  4.stl
    STL_files/SCS Lead Insulator.stl

They are still the published geometry, still in the model's own coordinate
frame, and they come back into any open document in one click --
build_lead_config.preview_rado_lead() imports them with no transform, into their
own group, as a disposable overlay that the lead designer panel has a button
for. So RADO's lead stops being a permanent conductor in the anatomy and becomes
what it actually is here: a reference placement to compare against.

It is worth having as exactly that. It is the only independently placed DRG lead
in the model -- RADO put it there by hand in SolidWorks -- and it is what the
foraminal corridor derived by measure_foramen.py is validated against
(`build_lead_config.py --name rado_drg_L3_match --compare-rado`). Deleting the
STLs, or the catalogue entry that reproduces it, would throw that away. Deleting
the five bodies from the document does not.

THE TRADEOFF, STATED PLAINLY
----------------------------
Against: NBF_RADO-SCS.FCStd stops being a faithful mirror of what RADO
published. Someone who opens it will not see the lead the paper's figures show,
and will have to know to press the button. That is a real cost and it is the
reason this is a separate, deliberate script rather than something that happens
quietly during some other job.
For: the document becomes usable as the mesh/solve input it is meant to be,
without a step where somebody has to remember to delete five bodies first -- and
the step people forget is the one that silently corrupts a solve rather than
failing it.

RUN IT FROM A RUNNING FREECAD, NOT freecadcmd
---------------------------------------------
This script SAVES the document, so it must not run headless. A .FCStd written
with no GUI is written without GuiDocument.xml, and every ShapeAppearance blob
in it -- all 245 tissue colours and transparencies -- is dropped silently
(verified on FreeCAD 26.3.0, git 48502). It has cost two recoveries already. So
main() refuses outright unless FreeCAD.GuiUp is true.

    # in FreeCAD's Python console, with NBF_RADO-SCS.FCStd open:
    import sys; sys.path.insert(0, "/home/mohamed/Projects/SCS-Modeling/src/freecad")
    import drop_rado_lead; drop_rado_lead.main()

    # dry run, no FreeCAD needed, reports the plan and touches nothing:
    python3 src/freecad/drop_rado_lead.py --dry-run

It is idempotent: run on a document that has already had them removed, it says
so and changes nothing.

Options come from DROP_RADO_LEAD_ARGS under a FreeCAD interpreter, for the
reason apply_labels.py documents: FreeCAD eats the command line itself. Output
goes to a report file as well as to stdout.
"""
import argparse
import os
import shlex
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
DEFAULT_DOC = os.path.join(REPO, "NBF_RADO-SCS.FCStd")
DEFAULT_REPORT = os.path.join(HERE, "drop_rado_lead_report.txt")

# RADO's own 4-contact DRG lead. Matched on Name, which FreeCAD makes immutable;
# Labels are apply_labels.py's and are not load-bearing.
RADO_LEAD = ("SCS_Lead_Electrode_1", "SCS_Lead_Electrode_2", "SCS_Lead_Electrode_3",
             "SCS_Lead_Electrode_4", "SCS_Lead_Insulator")

# The files those five bodies were imported from, which is why deleting them
# costs nothing. Checked to exist before anything is removed.
RADO_LEAD_STLS = ("SCS Lead Electrode 1.stl", "SCS Lead Electrode 2.stl",
                  "SCS Lead Electrode 3.stl", "SCS Lead Electrode 4.stl",
                  "SCS Lead Insulator.stl")


def script_args():
    """Arguments meant for this script; see apply_colors.script_args."""
    if os.path.basename(sys.argv[0]).lower().startswith("freecad"):
        return shlex.split(os.environ.get("DROP_RADO_LEAD_ARGS", ""))
    return sys.argv[1:]


def check_stls_present(stl_dir, say):
    """Refuse to delete the bodies unless the STLs they came from are on disk.

    The whole argument for removing them is that the geometry survives in
    STL_files/. If it does not, the argument does not hold and neither should
    the deletion.
    """
    missing = [fn for fn in RADO_LEAD_STLS
               if not os.path.exists(os.path.join(stl_dir, fn))]
    if missing:
        say("REFUSING: these bodies are only safe to delete because RADO's own")
        say("          STLs remain in %s, and %d of them are not there:",
            os.path.relpath(stl_dir, REPO), len(missing))
        for fn in missing:
            say("             %s", fn)
        return False
    say("source STLs: all %d present in %s -- the geometry survives this",
        len(RADO_LEAD_STLS), os.path.relpath(stl_dir, REPO))
    return True


def drop_rado_lead(doc, say):
    """Delete RADO's four contacts and its insulator from `doc`. Idempotent.

    A body is removed from its tissue group first. App::DocumentObjectGroup owns
    its children through a PropertyLinkList, and removing the object without
    removing the link is how a group ends up holding a dangling entry.
    """
    dropped, absent = [], []
    for name in RADO_LEAD:
        obj = doc.getObject(name)
        if obj is None:                      # already gone: a rerun, or a
            absent.append(name)              # document that never had it
            continue
        parent = obj.getParentGroup() if hasattr(obj, "getParentGroup") else None
        if parent is not None:
            parent.removeObject(obj)
        dropped.append("%s (%s)" % (name, obj.Label))
        doc.removeObject(name)
    say("deleted:   %d bodies -- %s", len(dropped), ", ".join(dropped) or "none")
    if absent:
        say("           already absent: %s", ", ".join(absent))
    return len(dropped)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--document", default=DEFAULT_DOC)
    ap.add_argument("--stl-dir", default=os.path.join(REPO, "STL_files"))
    ap.add_argument("--report", default=DEFAULT_REPORT)
    ap.add_argument("--dry-run", action="store_true",
                    help="report the plan; touches neither the document nor FreeCAD")
    args = ap.parse_args(script_args())

    out = []

    def say(fmt, *a):
        out.append(fmt % a if a else fmt)

    say("document:  %s", args.document)
    say("bodies:    %s", ", ".join(RADO_LEAD))
    say("")
    ok = check_stls_present(args.stl_dir, say)
    say("")

    rc = 0
    if args.dry_run:
        say("DRY RUN: would remove those %d bodies and save the document.",
            len(RADO_LEAD))
        say("They come back into any open document as a disposable overlay via")
        say("build_lead_config.preview_rado_lead(), or the lead designer panel's")
        say("\"Show RADO's lead\" button.")
    elif not ok:
        rc = 2
    else:
        try:
            import FreeCAD
        except ImportError:
            say("This needs to run inside FreeCAD. Use --dry-run here, or run")
            say("main() from FreeCAD's Python console with the document open.")
            rc = 2
            FreeCAD = None
        if rc == 0:
            if not FreeCAD.GuiUp:
                # The rule this whole pipeline is built around; see the docstring.
                say("REFUSING: FreeCAD.GuiUp is false. Saving a document without a")
                say("          GUI drops GuiDocument.xml and every ShapeAppearance")
                say("          blob in it -- all 245 tissue colours, silently.")
                say("          Run this from a running FreeCAD instead.")
                rc = 2
            else:
                doc = None
                for candidate in FreeCAD.listDocuments().values():
                    if os.path.abspath(getattr(candidate, "FileName", "")) == \
                            os.path.abspath(args.document):
                        doc = candidate
                        break
                if doc is None:
                    doc = FreeCAD.openDocument(args.document)
                    say("opened:    %s", doc.Name)
                else:
                    say("using the already-open document %s", doc.Name)
                before = len(doc.Objects)
                n = drop_rado_lead(doc, say)
                doc.recompute()
                if n:
                    doc.save()
                    say("saved:     %s, %d objects (was %d)",
                        args.document, len(doc.Objects), before)
                else:
                    say("nothing to do: the document already holds anatomy only.")

    say("")
    say("RESULT: %s", {0: "dry run OK" if args.dry_run else "OK",
                       2: "refused"}[rc])
    with open(args.report, "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    sys.stdout.write("\n".join(out) + "\n")
    return rc


def run_as_freecadcmd_script():
    """True when a FreeCAD interpreter was handed THIS file to run; see apply_labels.py."""
    if len(sys.argv) < 2 or not os.path.basename(sys.argv[0]).lower().startswith("freecad"):
        return False
    return os.path.abspath(sys.argv[1]) == os.path.abspath(__file__)


if __name__ == "__main__":
    sys.exit(main())
elif run_as_freecadcmd_script():
    main()
