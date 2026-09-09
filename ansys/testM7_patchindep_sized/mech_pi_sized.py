# Test M7: patch-independent tet mesh of a faceted STL body, with the sizing
# inputs that Mechanical actually requires.
#
# Trail that got here:
#   M2  headless Mechanical runs and imports STL, but mesh = 0 nodes
#   M3  produced NO output at all -> the -script engine is IronPython 2.7.4,
#       and py3-only syntax is a silent compile error
#   M4  body is a real solid (GeoBodySolid, 8835.8 mm^3) but faces=1, edges=0,
#       vertices=0 -> tessellated, no feature topology; default patch-CONFORMING
#       mesher reports "The mesh generation did not complete."
#   M5  Method=AllTriAllTet + Algorithm=MeshMethodAlgorithm.PatchIndependent is
#       accepted (AlgorithmType was the wrong enum -- that one is CMFD/MFD/
#       ProgramControlled/SCPIP), but meshing fails with "A mesh control or
#       method control does not have some required input defined."
#   M6  enumerated the method's 377 properties: the sizing knobs are named
#       MaximumElementSize / MinimumSizeLimit / DefeaturingTolerance /
#       CurvatureNormalAngle -- NOT ElementSize/MaxElementSize/MinElementSize.
#       MaximumElementSize was 0 [m] = the undefined required input.
#
# IronPython 2.7 only.
import os
import sys
import traceback


def log(m):
    print("[M7] " + str(m))
    sys.stdout.flush()


def messages(tag):
    try:
        msgs = ExtAPI.Application.Messages
        log("--- messages (" + tag + "): " + str(msgs.Count) + " ---")
        for i in range(msgs.Count):
            log("    [" + str(msgs[i].Severity) + "] " + str(msgs[i].DisplayString))
    except Exception:
        log("messages failed")


STL = os.path.expanduser("~/Documents/SCS-Modeling/NBF_RADO-SCS_STL_corrected/T8-10 - neuro_CSF-1.STL")
OUT = os.path.expanduser("~/Documents/SCS-Modeling/ansys_training/testM7_patchindep_sized")

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
log("body " + str(body.Name) + " volume " + str(body.Volume))

mesh = model.Mesh
sel = ExtAPI.SelectionManager.CreateSelectionInfo(
    Ansys.ACT.Interfaces.Common.SelectionTypeEnum.GeometryEntities)
sel.Ids = [body.GetGeoBody().Id]

meth = mesh.AddAutomaticMethod()
meth.Location = sel
meth.Method = E.MethodType.AllTriAllTet
meth.Algorithm = E.MeshMethodAlgorithm.PatchIndependent

# The names that actually exist (see M6). MaximumElementSize is the required one.
settings = (
    ("MaximumElementSize", "2 [mm]"),
    ("MinimumSizeLimit", "0.4 [mm]"),
    ("DefeaturingTolerance", "0.05 [mm]"),
    ("CurvatureNormalAngle", "18 [deg]"),
    ("FeatureAngle", "30 [deg]"),
)
for name, val in settings:
    try:
        setattr(meth, name, Quantity(val))
        log("meth." + name + " = " + val)
    except Exception:
        log("meth." + name + " FAILED: " + traceback.format_exc().splitlines()[-1])

for name, val in (("MeshBasedDefeaturing", True), ("RefineSurfaceMesh", True)):
    try:
        setattr(meth, name, val)
        log("meth." + name + " = " + str(val))
    except Exception:
        log("meth." + name + " FAILED")

try:
    mesh.GenerateMesh()
    log("RESULT: nodes=" + str(mesh.Nodes) + " elements=" + str(mesh.Elements))
except Exception:
    log("GenerateMesh raised: " + traceback.format_exc().splitlines()[-1])
messages("after mesh")

if mesh.Nodes > 0:
    log("SUCCESS -- patch-independent tets mesh a faceted STL body")
    try:
        mesh.MeshMetric = E.MeshMetricType.ElementQuality
        log("element quality: min=" + str(mesh.MeshMetricMin)
            + " avg=" + str(mesh.MeshMetricAverage)
            + " max=" + str(mesh.MeshMetricMax))
    except Exception:
        log("mesh metric unavailable: " + traceback.format_exc().splitlines()[-1])
    try:
        if not os.path.isdir(OUT):
            os.makedirs(OUT)
        ExtAPI.DataModel.Project.SaveAs(os.path.join(OUT, "m7.mechdat"), True)
        log("saved " + os.path.join(OUT, "m7.mechdat"))
    except Exception:
        log("save failed: " + traceback.format_exc().splitlines()[-1])
else:
    log("STILL 0 NODES -- see messages above")
log("done")
