"""Export the solved field in the forms downstream stages need.

  solution.vtu      unstructured VTK -- opens in ParaView and in FreeCAD's FEM
                    post-processing pipeline (FEM -> Post -> open result file)
  mesh_tagged.msh   gmsh 2.2 with ONE PHYSICAL VOLUME PER TISSUE, so the same
                    mesh can be handed to ElmerGrid/ElmerSolver (or FreeCAD's
                    Elmer StatCurrentSolver) once Elmer is installed, with the
                    material regions already separated
  field_grid.npz    structured phi grid for the NEURON stage (see field.py)
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
import field as F

GRID_SPACING = 0.25   # mm, the coarse end of neuron_plan chunk 3's 0.1-0.25 mm


def write_vtu(d, path):
    import meshio
    lab = d["label"]
    cells = [("tetra", d["tets"])]
    meshio.write_points_cells(
        path, d["nodes"], cells,
        point_data={"V_volts": d["V"].astype(np.float64),
                    **{"phi_c%d_ohm" % (i + 1): d["phi"][i].astype(np.float64)
                       for i in range(8)}},
        cell_data={"sigma_S_per_m": [d["sigma"]],
                   "tissue_id": [lab.astype(np.int32)]})


def write_tagged_msh(d, path):
    nodes, tets, lab = d["nodes"], d["tets"], d["label"]
    order = [str(x) for x in d["order"]]
    used = sorted(set(int(x) for x in np.unique(lab)))
    with open(path, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        f.write("$PhysicalNames\n%d\n" % len(used))
        for k, u in enumerate(used, start=1):
            f.write('3 %d "%s"\n' % (k, order[u] if u >= 0 else "background"))
        f.write("$EndPhysicalNames\n")
        f.write("$Nodes\n%d\n" % len(nodes))
        np.savetxt(f, np.column_stack([np.arange(1, len(nodes) + 1), nodes]),
                   fmt="%d %.9g %.9g %.9g")
        f.write("$EndNodes\n$Elements\n%d\n" % len(tets))
        tag = np.searchsorted(used, lab) + 1
        arr = np.column_stack([np.arange(1, len(tets) + 1),
                               np.full(len(tets), 4), np.full(len(tets), 2),
                               tag, tag, tets + 1])
        np.savetxt(f, arr, fmt="%d")
        f.write("$EndElements\n")


def main():
    sol = os.path.join(C.OUT, "solution.npz")
    d = np.load(sol, allow_pickle=True)
    write_vtu(d, os.path.join(C.OUT, "solution.vtu"))
    print("wrote solution.vtu")
    write_tagged_msh(d, os.path.join(C.OUT, "mesh_tagged.msh"))
    print("wrote mesh_tagged.msh")

    # structured grid over cord + surrounding sac, spanning the lead generously
    # Cover the whole canal cross-section (not just the cord) over the lead's
    # rostrocaudal span plus margin, clipped to the model's own z extent
    # (RADO stops at z = 164.9 mm). Points outside the tissue union come back
    # NaN, which is the honest answer -- there is no tissue there to sample.
    lo = np.array([45.0, 71.5, 104.0])
    hi = np.array([68.0, 90.5, 164.5])
    shape = tuple(int(np.floor((hi[i] - lo[i]) / GRID_SPACING)) + 1 for i in range(3))
    cov = F.structured_export(sol, os.path.join(C.OUT, "field_grid.npz"),
                              lo, GRID_SPACING, shape)
    print("wrote field_grid.npz  shape=%s  spacing=%.2f mm  in-mesh coverage %.1f %%"
          % (shape, GRID_SPACING, 100 * cov))

    # sanity check demanded by neuron_plan chunk 3: phi must fall off with r
    g = np.load(os.path.join(C.OUT, "field_grid.npz"), allow_pickle=True)
    phi = g["phi"][C.SOURCE_CONTACT - 1]
    ax = [lo[i] + GRID_SPACING * np.arange(shape[i]) for i in range(3)]
    X, Y, Z = np.meshgrid(*ax, indexing="ij")
    c = np.array([55.60, 87.49, 132.83])          # contact-3 centroid
    r = np.sqrt((X - c[0]) ** 2 + (Y - c[1]) ** 2 + (Z - c[2]) ** 2)
    m = np.isfinite(phi) & (r > 1.5) & (r < 6.0)
    a = np.polyfit(np.log(r[m]), np.log(np.abs(phi[m] - np.nanmedian(phi))), 1)[0]
    print("  phi(r) around contact 3 falls as r^%.2f over 1.5-6 mm "
          "(a point source in a uniform medium would give -1.00; the layered "
          "anatomy is expected to deviate)" % a)


if __name__ == "__main__":
    main()
