# Test M2: why did GenerateMesh() produce 0 nodes for an imported STL?
# Inspect what the STL actually became (solid vs surface vs faceted), read the
# Mechanical message log, then try an explicit patch-independent tet method.
import os, sys, traceback

def log(m):
    print("[M2] %s" % m); sys.stdout.flush()

def dump_messages(tag):
    try:
        msgs = ExtAPI.DataModel.Project.Messages
        log("--- messages (%s): %d ---" % (tag, msgs.Count))
        for i in range(msgs.Count):
            m = msgs[i]
            log("    [%s] %s | %s" % (m.Severity, m.DisplayString, getattr(m, "Location", "")))
    except Exception:
        log("messages read failed:\n" + traceback.format_exc())

STL = os.path.expanduser("~/Documents/SCS-Modeling/NBF_RADO-SCS_STL_corrected/T8-10 - neuro_CSF-1.STL")
OUT = os.path.expanduser("~/Documents/SCS-Modeling/ansys_training/testM2_mechanical_mesh")

model = ExtAPI.DataModel.Project.Model
log("ProductVersion %s" % ExtAPI.DataModel.Project.ProductVersion)

# ---- import -----------------------------------------------------------------
try:
    gi = model.GeometryImportGroup.AddGeometryImport()
    fmt = Ansys.Mechanical.DataModel.Enums.GeometryImportPreference.Format.Automatic
    prefs = Ansys.ACT.Mechanical.Utilities.GeometryImportPreferences()
    prefs.ProcessNamedSelections = True
    prefs.ProcessSolids = True
    gi.Import(STL, fmt, prefs)
    log("import ok")
except Exception:
    log("import failed:\n" + traceback.format_exc())

dump_messages("after import")

# ---- what did we actually get? ---------------------------------------------
bodies = []
try:
    geo = model.Geometry
    log("Geometry.Children = %d" % geo.Children.Count)
    for part in geo.Children:
        log("  Part '%s' children=%d" % (part.Name, part.Children.Count))
        for b in part.Children:
            bodies.append(b)
            info = []
            for attr in ("Name", "Suppressed", "Dimension", "StiffnessBehavior", "Volume", "Visible"):
                try: info.append("%s=%s" % (attr, getattr(b, attr)))
                except Exception: info.append("%s=<n/a>" % attr)
            log("    Body: " + ", ".join(str(x) for x in info))
            try:
                gb = b.GetGeoBody()
                log("      GeoBody type=%s cells=%s faces=%s" % (
                    type(gb).__name__,
                    len(gb.Cells) if hasattr(gb, "Cells") else "?",
                    len(gb.Faces) if hasattr(gb, "Faces") else "?"))
            except Exception:
                log("      GetGeoBody failed: " + traceback.format_exc().splitlines()[-1])
except Exception:
    log("geometry walk failed:\n" + traceback.format_exc())

# ---- mesh attempt 1: default ------------------------------------------------
mesh = model.Mesh
try:
    mesh.ElementSize = Quantity("1 [mm]")
    mesh.GenerateMesh()
    log("default mesh -> nodes=%s elements=%s" % (mesh.Nodes, mesh.Elements))
except Exception:
    log("default mesh failed:\n" + traceback.format_exc())
dump_messages("after default mesh")

# ---- mesh attempt 2: explicit patch-independent tetrahedrons ----------------
try:
    sel = ExtAPI.SelectionManager.CreateSelectionInfo(
        Ansys.ACT.Interfaces.Common.SelectionTypeEnum.GeometryEntities)
    ids = []
    for b in bodies:
        try: ids.append(b.GetGeoBody().Id)
        except Exception: pass
    log("body ids for method scoping: %s" % ids)
    sel.Ids = ids
    meth = mesh.AddAutomaticMethod()
    meth.Location = sel
    E = Ansys.Mechanical.DataModel.Enums
    meth.Method = E.MethodType.Tetrahedrons
    meth.Algorithm = E.AlgorithmType.PatchIndependent
    log("patch-independent tet method added; regenerating...")
    mesh.GenerateMesh()
    log("patch-independent mesh -> nodes=%s elements=%s" % (mesh.Nodes, mesh.Elements))
except Exception:
    log("patch-independent mesh failed:\n" + traceback.format_exc())
dump_messages("after patch-independent mesh")

# ---- save -------------------------------------------------------------------
try:
    if not os.path.isdir(OUT): os.makedirs(OUT)
    ExtAPI.DataModel.Project.SaveAs(os.path.join(OUT, "m2.mechdat"), True)
    log("saved")
except Exception:
    log("save failed:\n" + traceback.format_exc())
log("done")
