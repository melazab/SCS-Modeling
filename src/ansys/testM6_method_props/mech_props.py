# Test M6: enumerate the real property names on a PatchIndependent AutomaticMethod
# (and on the global Mesh) instead of guessing them.
#
# M5 got Method=AllTriAllTet + Algorithm=PatchIndependent accepted, but meshing
# failed with:
#     [Error] A mesh control or method control does not have some required
#             input defined. Please define.
# and ElementSize / MaxElementSize / MinElementSize were all rejected with
# AttributeError on the AutomaticMethod object. Mechanical ACT objects expose a
# Properties collection whose entries carry the true names, so dump that.
#
# IronPython 2.7 only.
import os
import sys
import traceback


def log(m):
    print("[M6] " + str(m))
    sys.stdout.flush()


def dump_props(obj, label):
    log("=== properties of " + label + " ===")
    for coll_name in ("Properties", "VisibleProperties"):
        try:
            coll = getattr(obj, coll_name)
        except Exception:
            log("  (no " + coll_name + ")")
            continue
        try:
            n = coll.Count
        except Exception:
            try:
                n = len(coll)
            except Exception:
                log("  (" + coll_name + " not countable)")
                continue
        log("  -- " + coll_name + ": " + str(n) + " --")
        for i in range(n):
            try:
                p = coll[i]
                name = str(getattr(p, "Name", "?"))
                try:
                    val = str(p.Value)
                except Exception:
                    val = "<unreadable>"
                ro = ""
                try:
                    if p.ReadOnly:
                        ro = " [read-only]"
                except Exception:
                    pass
                log("     " + name + " = " + val + ro)
            except Exception:
                log("     <property " + str(i) + " failed>")


def messages(tag):
    try:
        msgs = ExtAPI.Application.Messages
        log("--- messages (" + tag + "): " + str(msgs.Count) + " ---")
        for i in range(msgs.Count):
            log("    [" + str(msgs[i].Severity) + "] " + str(msgs[i].DisplayString))
    except Exception:
        log("messages failed")


STL = os.path.expanduser("~/Documents/SCS-Modeling/NBF_RADO-SCS_STL_corrected/T8-10 - neuro_CSF-1.STL")
OUT = os.path.expanduser("~/Documents/SCS-Modeling/ansys_training/testM6_method_props")

E = Ansys.Mechanical.DataModel.Enums
model = ExtAPI.DataModel.Project.Model
gi = model.GeometryImportGroup.AddGeometryImport()
prefs = Ansys.ACT.Mechanical.Utilities.GeometryImportPreferences()
prefs.ProcessNamedSelections = True
prefs.ProcessSolids = True
gi.Import(STL, E.GeometryImportPreference.Format.Automatic, prefs)

body = None
for part in model.Geometry.Children:
    for b in part.Children:
        body = b
log("body " + str(body.Name))

mesh = model.Mesh
dump_props(mesh, "Mesh (global)")

sel = ExtAPI.SelectionManager.CreateSelectionInfo(
    Ansys.ACT.Interfaces.Common.SelectionTypeEnum.GeometryEntities)
sel.Ids = [body.GetGeoBody().Id]
meth = mesh.AddAutomaticMethod()
meth.Location = sel
meth.Method = E.MethodType.AllTriAllTet
meth.Algorithm = E.MeshMethodAlgorithm.PatchIndependent
dump_props(meth, "AutomaticMethod (AllTriAllTet + PatchIndependent)")

log("=== settable-looking attributes on the method ===")
names = []
for a in dir(meth):
    if a.startswith("_"):
        continue
    low = a.lower()
    if ("size" in low or "angle" in low or "defeat" in low or "curv" in low
            or "prox" in low or "growth" in low or "refine" in low or "tol" in low):
        names.append(a)
log(str(names))
for a in names:
    try:
        log("  meth." + a + " = " + str(getattr(meth, a)))
    except Exception:
        log("  meth." + a + " <unreadable>")

log("done")
