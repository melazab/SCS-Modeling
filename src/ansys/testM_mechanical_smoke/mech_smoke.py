# Mechanical (standalone) batch smoke test. Run via:
#   /usr/local/ansys_inc/v251/aisol/.workbench -DSApplet -AppModeMech -b -script mech_smoke.py
# Goal: learn whether headless Mechanical works on Pioneer at all, whether it can
# import an STL directly, and whether it can mesh it. Every step is wrapped so
# the log shows exactly where it dies.
import os, sys, traceback

def log(msg):
    print("[SMOKE] " + str(msg))
    sys.stdout.flush()

STL = os.path.expanduser("~/Documents/SCS-Modeling/NBF_RADO-SCS_STL_corrected/T8-10 - neuro_CSF-1.STL")
OUT = os.path.expanduser("~/Documents/SCS-Modeling/ansys_training/testM_mechanical_smoke")

try:
    log("ProductVersion: %s" % ExtAPI.DataModel.Project.ProductVersion)
except Exception:
    log("ProductVersion failed:\n" + traceback.format_exc())

try:
    model = ExtAPI.DataModel.Project.Model
    log("Model: %s" % model.Name)
except Exception:
    log("Model access failed:\n" + traceback.format_exc())

try:
    log("Importing STL: %s (exists=%s)" % (STL, os.path.exists(STL)))
    gi_group = model.GeometryImportGroup
    gi = gi_group.AddGeometryImport()
    fmt = Ansys.Mechanical.DataModel.Enums.GeometryImportPreference.Format.Automatic
    prefs = Ansys.ACT.Mechanical.Utilities.GeometryImportPreferences()
    prefs.ProcessNamedSelections = True
    gi.Import(STL, fmt, prefs)
    geo = model.Geometry
    log("Geometry children: %d" % geo.Children.Count)
    for child in geo.Children:
        log("  part: %s (children %d)" % (child.Name, child.Children.Count))
except Exception:
    log("STL import failed:\n" + traceback.format_exc())

try:
    mesh = model.Mesh
    mesh.ElementSize = Quantity("0.5 [mm]")
    log("Generating mesh...")
    mesh.GenerateMesh()
    log("Mesh: %d nodes, %d elements" % (mesh.Nodes, mesh.Elements))
except Exception:
    log("Meshing failed:\n" + traceback.format_exc())

try:
    if not os.path.isdir(OUT):
        os.makedirs(OUT)
    path = os.path.join(OUT, "smoke.mechdat")
    ExtAPI.DataModel.Project.Save(path)
    log("Saved %s" % path)
except Exception:
    log("Save failed:\n" + traceback.format_exc())

log("done")
