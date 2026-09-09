# Test M3: stop guessing the Mechanical scripting API — enumerate it, and find
# out WHY an imported STL solid meshes to 0 nodes.
#
# Known from M2: the STL imports as a real 3D body (Dimension=Three_D,
# Volume=8835.8 mm^3, faces=1, not suppressed) but GenerateMesh() yields
# 0 nodes / 0 elements. Also: ExtAPI.DataModel.Project.Messages does NOT exist,
# and Ansys.Mechanical.DataModel.Enums.MethodType has no 'Tetrahedrons'.
import os, sys, traceback

def log(m):
    print("[M3] %s" % m); sys.stdout.flush()

def try_(label, fn):
    try:
        r = fn()
        log("%s -> %s" % (label, r))
        return r
    except Exception:
        log("%s FAILED: %s" % (label, traceback.format_exc().splitlines()[-1]))
        return None

STL = os.path.expanduser("~/Documents/SCS-Modeling/NBF_RADO-SCS_STL_corrected/T8-10 - neuro_CSF-1.STL")
OUT = os.path.expanduser("~/Documents/SCS-Modeling/ansys_training/testM3_mechanical_api")

# ---- 1. where do messages actually live? -----------------------------------
log("=== message API discovery ===")
for expr, fn in [
    ("ExtAPI.Application.Messages", lambda: ExtAPI.Application.Messages),
    ("ExtAPI.DataModel.Project.Model.Messages", lambda: ExtAPI.DataModel.Project.Model.Messages),
    ("ExtAPI.Application.MessageManager", lambda: ExtAPI.Application.MessageManager),
]:
    try_(expr, fn)

def messages(tag):
    try:
        msgs = ExtAPI.Application.Messages
        log("--- messages (%s): %d ---" % (tag, msgs.Count))
        for i in range(msgs.Count):
            m = msgs[i]
            log("    [%s] %s" % (m.Severity, m.DisplayString))
    except Exception:
        log("(messages unavailable: %s)" % traceback.format_exc().splitlines()[-1])

# ---- 2. enumerate the enums we need ----------------------------------------
log("=== enum discovery ===")
E = Ansys.Mechanical.DataModel.Enums
for enum_name in ("MethodType", "AlgorithmType", "MeshMethodAlgorithm", "ElementOrder",
                  "GeometryType", "MeshControlGroupType", "MethodMeshType"):
    try:
        e = getattr(E, enum_name)
        vals = [v for v in dir(e) if not v.startswith("_")]
        log("%s: %s" % (enum_name, vals))
    except Exception:
        log("%s: <not present>" % enum_name)

# ---- 3. import and inspect --------------------------------------------------
log("=== import ===")
model = ExtAPI.DataModel.Project.Model
gi = model.GeometryImportGroup.AddGeometryImport()
fmt = Ansys.Mechanical.DataModel.Enums.GeometryImportPreference.Format.Automatic
prefs = Ansys.ACT.Mechanical.Utilities.GeometryImportPreferences()
prefs.ProcessNamedSelections = True
prefs.ProcessSolids = True
gi.Import(STL, fmt, prefs)
messages("after import")

bodies = []
for part in model.Geometry.Children:
    for b in part.Children:
        bodies.append(b)
log("bodies: %d" % len(bodies))
b0 = bodies[0]
log("--- body attribute dump ---")
for attr in sorted(a for a in dir(b0) if not a.startswith("_")):
    if attr in ("Children", "Parent", "Properties", "VisibleProperties"):
        continue
    try:
        v = getattr(b0, attr)
        if callable(v):
            continue
        log("    %s = %s" % (attr, v))
    except Exception:
        pass

# is it a faceted/mesh body rather than a BREP solid?
gb = try_("body.GetGeoBody()", lambda: b0.GetGeoBody())
if gb is not None:
    for attr in ("Id", "Volume", "Area", "IsSolid", "BodyType", "Faces", "Edges", "Vertices"):
        try_("  geoBody.%s" % attr, lambda a=attr: getattr(gb, a))

# ---- 4. mesh, and read the messages this time -------------------------------
log("=== mesh attempts ===")
mesh = model.Mesh
try_("mesh.ElementSize=2mm", lambda: setattr(mesh, "ElementSize", Quantity("2 [mm]")))
try_("mesh.GenerateMesh()", lambda: mesh.GenerateMesh())
log("after default: nodes=%s elements=%s" % (mesh.Nodes, mesh.Elements))
messages("after default mesh")

# dump whatever mesh-state properties exist
for attr in ("MeshMetric", "CurrentConfiguration", "Worksheet", "ElementOrder",
             "PhysicsPreference", "ElementSize", "Resolution", "UseAdaptiveSizing"):
    try_("mesh.%s" % attr, lambda a=attr: getattr(mesh, a))

# patch-independent using whatever the enum is really called
try:
    sel = ExtAPI.SelectionManager.CreateSelectionInfo(
        Ansys.ACT.Interfaces.Common.SelectionTypeEnum.GeometryEntities)
    sel.Ids = [gb.Id]
    meth = mesh.AddAutomaticMethod()
    meth.Location = sel
    log("method defaults: Method=%s" % meth.Method)
    log("method Method options: %s" % [v for v in dir(type(meth.Method)) if not v.startswith("_")])
    messages("after adding method")
except Exception:
    log("method scoping failed:\n" + traceback.format_exc())

try:
    if not os.path.isdir(OUT): os.makedirs(OUT)
    ExtAPI.DataModel.Project.SaveAs(os.path.join(OUT, "m3.mechdat"), True)
    log("saved")
except Exception:
    log("save failed: %s" % traceback.format_exc().splitlines()[-1])
log("done")
