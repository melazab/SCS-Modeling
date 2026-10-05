"""Rescale the Elmer StatCurrentSolver result from metres back to millimetres.

`crosscheck_elmer.py` runs Elmer with `Coordinate Scaling = 1.0e-3` in its
.sif (needed because the mesh was built in mm, matching the RADO STL files,
but Elmer's own convention is metres). That means the point coordinates
stored in the output VTU,

    fem/out/elmer/results/case_t0001.vtu

are in **metres**, while every other artifact in this repo (STLs, the P1
solver's fem/out/solution.vtu, the FreeCAD document) is in **millimetres**.
Loading case_t0001.vtu directly into FreeCAD's Post pipeline puts the field
~1000x too small and mis-positioned relative to the anatomy.

This script writes a millimetre-scale copy alongside it:

    fem/out/elmer/results/case_t0001_mm.vtu

by reading the VTU with vtk's own XML reader, multiplying every point
coordinate by 1000, and writing it back out with vtk's own XML writer --
point data (`potential`, `volume current`) and cell data (`GeometryIds`) pass
through untouched, only the geometry moves. No physics, mesh connectivity,
or field values are altered.

Run whenever a new Elmer solve produces a new case_t0001.vtu:

    python3 fem/scripts/prepare_elmer_result_mm.py

`fem/scripts/view_elmer_result.FCMacro` calls this automatically (via
`prepare()`) if the mm file is missing or older than the source, so this
does not normally need to be run by hand.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C

SRC = os.path.join(C.OUT, "elmer", "results", "case_t0001.vtu")
DST = os.path.join(C.OUT, "elmer", "results", "case_t0001_mm.vtu")

SCALE = 1000.0  # metres -> millimetres


def prepare(src=SRC, dst=DST, force=False):
    """Write the mm-scaled copy at `dst` if missing or stale. Returns dst."""
    if not force and os.path.exists(dst) and os.path.getmtime(dst) >= os.path.getmtime(src):
        print("up to date:", dst)
        return dst

    import vtk

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(src)
    reader.Update()
    grid = reader.GetOutput()

    n_pts = grid.GetNumberOfPoints()
    n_cells = grid.GetNumberOfCells()
    if n_pts == 0:
        raise RuntimeError("read 0 points from %s -- is it a valid VTU?" % src)

    bounds_m = grid.GetBounds()

    pts = grid.GetPoints()
    from vtk.util import numpy_support as vnp
    arr = vnp.vtk_to_numpy(pts.GetData()).copy()
    arr *= SCALE
    new_arr = vnp.numpy_to_vtk(arr, deep=True)
    new_arr.SetName(pts.GetData().GetName())
    pts.SetData(new_arr)
    grid.Modified()

    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(dst)
    writer.SetInputData(grid)
    writer.SetDataModeToAppended()
    writer.EncodeAppendedDataOff()  # raw binary, matches the source file's encoding
    writer.Write()

    print("read  %d points, %d cells from %s" % (n_pts, n_cells, src))
    print("bounds (m):  X%s Y%s Z%s" % (bounds_m[0:2], bounds_m[2:4], bounds_m[4:6]))
    bounds_mm = grid.GetBounds()
    print("bounds (mm): X%s Y%s Z%s" % (bounds_mm[0:2], bounds_mm[2:4], bounds_mm[4:6]))
    print("wrote", dst)
    return dst


if __name__ == "__main__":
    prepare(force=True)
