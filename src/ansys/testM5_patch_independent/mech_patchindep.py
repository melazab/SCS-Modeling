# Test M5: mesh a faceted STL body with a PATCH INDEPENDENT tetrahedral method.
#
# Why: M4 showed the imported STL is a real solid (GeoBodySolid, 8835.8 mm^3)
# but has faces=1, edges=0, vertices=0 -- a tessellated body with no feature
# topology. The default patch-CONFORMING mesher must respect every face patch,
# and on that geometry it just gives up:
#     [Error] The mesh generation did not complete.
# Patch independent tets ignore the surface patch structure and mesh the
# enclosed volume, which is the standard remedy for faceted/dirty geometry.
#
# M4 also showed AlgorithmType is the WRONG enum (its values are CMFD/MFD/
# ProgramControlled/SCPIP -- optimisation algorithms). The mesh method algorithm
# lives on MeshMethodAlgorithm.
#
# IronPython 2.7 only -- no py3 syntax (see the M3 silent-compile-error trap).
import os
import sys
import traceback


def log(m):
    print("[M5] " + str(m))
    sys.stdout.flush()


def enum_values(e):
    skip = ("CompareTo", "Equals", "Format", "GetHashCode", "GetName", "GetNames",
            "GetType", "GetTypeCode", "GetUnderlyingType", "GetValues", "HasFlag",
            "IsDefined", "Parse", "ToObject", "ToString", "value__",
            "MemberwiseClone", "ReferenceEquals")
    out = []
    for v in dir(e):
        if v.startswith("_") or v in skip or v.startswith("To"):
            continue
        out.append(v)
    return out


def messages(tag):
    try:
        msgs = ExtAPI.Application.Messages
        log("--- messages (" + tag + "): " + str(msgs.Count) + " ---")
        for i in range(msgs.Count):
            log("    [" + str(msgs[i].Severity) + "] " + str(msgs[i].DisplayString))
    except Exception:
        log("messages failed: " + traceback.format_exc().splitlines()[-1])


STL = os.path.expanduser("~/Documents/SCS-Modeling/NBF_RADO-SCS_STL_corrected/T8-10 - neuro_CSF-1.STL")
OUT = os.path.expanduser("~/Documents/SCS-Modeling/ansys_training/testM5_patch_independent")

E = Ansys.Mechanical.DataModel.Enums
for name in ("MeshMethodAlgorithm", "MethodType", "ElementOrder", "MeshControlGroupType"):
    try:
        log(name + ": " + str(enum_values(getattr(E, name))))
    except Exception:
        log(name + ": <not present>")

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
log("body: " + str(body.Name) + " volume=" + str(body.Volume))

mesh = model.Mesh
sel = ExtAPI.SelectionManager.CreateSelectionInfo(
    Ansys.ACT.Interfaces.Common.SelectionTypeEnum.GeometryEntities)
sel.Ids = [body.GetGeoBody().Id]

meth = mesh.AddAutomaticMethod()
meth.Location = sel
meth.Method = E.MethodType.AllTriAllTet
log("Method = " + str(meth.Method))

# The whole point of this test:
try:
    meth.Algorithm = E.MeshMethodAlgorithm.PatchIndependent
    log("Algorithm = " + str(meth.Algorithm))
except Exception:
    log("setting PatchIndependent failed: " + traceback.format_exc().splitlines()[-1])
    log("available on meth: " + str([a for a in dir(meth) if "lgorith" in a or "atch" in a]))

# Patch independent sizing controls. The CSF body is ~30 mm across and the
# thinnest structures elsewhere in this model are ~0.4 mm, so keep the max
# element modest and let it refine.
for attr, val in (("ElementSize", "2 [mm]"),
                  ("MaxElementSize", "2 [mm]"),
                  ("MinElementSize", "0.3 [mm]"),
                  ("FeatureAngle", "30 [deg]")):
    try:
        setattr(meth, attr, Quantity(val))
        log("meth." + attr + " = " + val)
    except Exception:
        log("meth." + attr + " not settable: " + traceback.format_exc().splitlines()[-1])

for attr, val in (("Defeaturing", True), ("MeshBasedDefeaturing", True)):
    try:
        setattr(meth, attr, val)
        log("meth." + attr + " = " + str(val))
    except Exception:
        pass

try:
    mesh.GenerateMesh()
    log("RESULT patch-independent: nodes=" + str(mesh.Nodes)
        + " elements=" + str(mesh.Elements))
except Exception:
    log("GenerateMesh raised: " + traceback.format_exc().splitlines()[-1])
messages("after patch-independent mesh")

# If it worked, write the mesh out so MAPDL can consume it.
if mesh.Nodes > 0:
    try:
        if not os.path.isdir(OUT):
            os.makedirs(OUT)
        ExtAPI.DataModel.Project.SaveAs(os.path.join(OUT, "m5.mechdat"), True)
        log("saved project with mesh")
    except Exception:
        log("save failed: " + traceback.format_exc().splitlines()[-1])
log("done")
