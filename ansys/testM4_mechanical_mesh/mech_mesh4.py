# Test M4: why does an imported STL solid mesh to 0 nodes?
#
# IMPORTANT: Mechanical's -script engine is IronPython 2.7.4 (on Mono), NOT
# CPython 3. Any py3-only syntax is a compile error and the script produces
# ZERO output (silently), which is what happened to Test M3. Keep this file
# strictly Python 2.7 compatible: no f-strings, no print(..., flush=), no
# dict/set comprehension tricks, no `nonlocal`.
import os
import sys
import traceback


def log(m):
    print("[M4] " + str(m))
    sys.stdout.flush()


def enum_values(e):
    out = []
    skip = ("CompareTo", "Equals", "Format", "GetHashCode", "GetName", "GetNames",
            "GetType", "GetTypeCode", "GetUnderlyingType", "GetValues", "HasFlag",
            "IsDefined", "Parse", "ToObject", "ToString", "value__")
    for v in dir(e):
        if v.startswith("_") or v in skip:
            continue
        out.append(v)
    return out


def messages(tag):
    try:
        msgs = ExtAPI.Application.Messages
        log("--- messages (" + tag + "): " + str(msgs.Count) + " ---")
        for i in range(msgs.Count):
            m = msgs[i]
            log("    [" + str(m.Severity) + "] " + str(m.DisplayString))
    except Exception:
        log("messages failed: " + traceback.format_exc().splitlines()[-1])


STL = os.path.expanduser("~/Documents/SCS-Modeling/NBF_RADO-SCS_STL_corrected/T8-10 - neuro_CSF-1.STL")
OUT = os.path.expanduser("~/Documents/SCS-Modeling/ansys_training/testM4_mechanical_mesh")

E = Ansys.Mechanical.DataModel.Enums
log("MethodType: " + str(enum_values(E.MethodType)))
log("AlgorithmType: " + str(enum_values(E.AlgorithmType)))

model = ExtAPI.DataModel.Project.Model

# ---- import -----------------------------------------------------------------
gi = model.GeometryImportGroup.AddGeometryImport()
fmt = E.GeometryImportPreference.Format.Automatic
prefs = Ansys.ACT.Mechanical.Utilities.GeometryImportPreferences()
prefs.ProcessNamedSelections = True
prefs.ProcessSolids = True
gi.Import(STL, fmt, prefs)
log("imported")
messages("after import")

bodies = []
for part in model.Geometry.Children:
    for b in part.Children:
        bodies.append(b)
log("body count: " + str(len(bodies)))
b0 = bodies[0]

for attr in ("Name", "Suppressed", "Dimension", "Volume", "Material",
             "StiffnessBehavior", "GeometryType", "MeshMetric"):
    try:
        log("body." + attr + " = " + str(getattr(b0, attr)))
    except Exception:
        log("body." + attr + " unavailable")

# Is it a faceted (tessellated) body rather than a BREP solid? That is the
# prime suspect for "meshes to 0 nodes".
try:
    gb = b0.GetGeoBody()
    log("geoBody: " + str(gb))
    for attr in ("Id", "Volume", "Area", "BodyType", "IsSolid"):
        try:
            log("  geoBody." + attr + " = " + str(getattr(gb, attr)))
        except Exception:
            log("  geoBody." + attr + " unavailable")
    try:
        log("  geoBody faces=" + str(gb.Faces.Count) + " edges=" + str(gb.Edges.Count)
            + " vertices=" + str(gb.Vertices.Count))
    except Exception:
        log("  face/edge/vertex counts unavailable: "
            + traceback.format_exc().splitlines()[-1])
except Exception:
    log("GetGeoBody failed: " + traceback.format_exc().splitlines()[-1])

mesh = model.Mesh

# ---- attempt 1: default automatic mesh --------------------------------------
try:
    mesh.ElementSize = Quantity("2 [mm]")
    mesh.GenerateMesh()
    log("attempt1 default: nodes=" + str(mesh.Nodes) + " elements=" + str(mesh.Elements))
except Exception:
    log("attempt1 failed: " + traceback.format_exc().splitlines()[-1])
messages("after attempt1")

# ---- attempt 2: explicit tet method, patch independent ----------------------
try:
    sel = ExtAPI.SelectionManager.CreateSelectionInfo(
        Ansys.ACT.Interfaces.Common.SelectionTypeEnum.GeometryEntities)
    sel.Ids = [b0.GetGeoBody().Id]
    meth = mesh.AddAutomaticMethod()
    meth.Location = sel
    log("method default Method=" + str(meth.Method))
    meth.Method = E.MethodType.AllTriAllTet          # 'Tetrahedrons' does NOT exist
    log("set Method=AllTriAllTet")
    try:
        meth.Algorithm = E.AlgorithmType.PatchIndependent
        log("set Algorithm=PatchIndependent")
    except Exception:
        log("PatchIndependent unavailable: " + traceback.format_exc().splitlines()[-1])
    mesh.GenerateMesh()
    log("attempt2 patch-independent: nodes=" + str(mesh.Nodes)
        + " elements=" + str(mesh.Elements))
except Exception:
    log("attempt2 failed: " + traceback.format_exc().splitlines()[-1])
messages("after attempt2")

try:
    if not os.path.isdir(OUT):
        os.makedirs(OUT)
    ExtAPI.DataModel.Project.SaveAs(os.path.join(OUT, "m4.mechdat"), True)
    log("saved")
except Exception:
    log("save failed: " + traceback.format_exc().splitlines()[-1])
log("done")
