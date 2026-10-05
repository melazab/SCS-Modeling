"""Build, classify and preview a live lead mesh in a scratch run directory.

CLI: mesh_preview.py <run_dir> <params.json>
Progress 0-75 is the volume mesh, 75-85 classification, 85-92 the dura check,
and 92-100 tissue-surface export. Face sorting is partitioned on disk for
bounded memory. A mesh_report and artifact manifests mark completed previews.
An interrupted preview can reuse a completed matching volume mesh.
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_mesh
import config as C
from assign_and_solve import classify_tets, dura_leak_report


from mesh_faces import face_batches
import artifacts


def outer_surface(nodes, tets, lab):
    """Every triangular face owned by exactly one tet (the mesh's outer
    boundary, tissue-labelled by whichever tet owns it), for a lightweight
    preview -- same face-counting idiom assign_and_solve.py already uses.

    SUPERSEDED for the preview by tissue_surfaces() below, which also emits
    the faces BETWEEN tissues. Kept because it is the honest answer to
    "what is the outer hull of this mesh", which is a different question.
    """
    faces, ids = [], []
    for f, owner, paired, lone in face_batches(tets):
        faces.append(f[lone])
        ids.append(lab[owner[lone]])
    return np.concatenate(faces), np.concatenate(ids)


def tissue_surfaces(nodes, tets, lab):
    """Every face that BOUNDS A TISSUE, labelled with the tissue it bounds.

    WHY THIS REPLACED outer_surface() FOR THE PREVIEW (2026-09-15)

    The preview used to be the mesh's outer HULL only -- the faces owned by
    exactly one tet. Mohamed, looking at one: "this mesh just looks like a
    blob to me". It genuinely is one, and no amount of extra resolution
    fixes it: the outer hull of a box mesh of the spinal canal is a closed
    envelope, so the only thing it can ever show is the outside of the
    outermost compartment. That mattered less when the preview was one fused
    object; it is fatal now that the point is one toggleable object per
    tissue, because every tissue except the outermost had almost nothing on
    that hull to draw. Measured on the field-of-record mesh: 39 916 hull
    facets against 348 871 faces that actually separate one tissue from
    another.

    So this returns both:
      * a face on the domain boundary, once, labelled with its owner's tissue;
      * a face between two tets of DIFFERENT tissues, TWICE -- once for each
        side -- so each tissue gets its own complete shell and hiding the
        epidural space in the tree reveals the dura underneath it, exactly
        the way hiding "04 Epidural space" reveals "05 Dura and meninges" in
        the anatomy;
      * nothing for a face between two tets of the same tissue (that is the
        interior, and drawing it would be both useless and enormous).

    Cost, measured: 39 916 -> 737 658 facets on the field-of-record mesh,
    against the 1 771 980 triangles RADO's own 240 STL bodies already render
    in this document. Returns (tris (m,3), tissue_id (m,)).
    """
    tris, ids = [], []
    for f, owner, paired, lone in face_batches(tets):
        tris.append(f[lone])
        ids.append(lab[owner[lone]])
        la, lb = lab[owner[paired]], lab[owner[paired + 1]]
        cut = la != lb
        tris.extend((f[paired[cut]], f[paired[cut] + 1]))
        ids.extend((la[cut], lb[cut]))
    return np.concatenate(tris), np.concatenate(ids)


def write_preview_vtu(path, nodes, tris, tissue_id):
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk

    pts = vtk.vtkPoints()
    pts.SetData(numpy_to_vtk(np.ascontiguousarray(nodes, dtype=np.float64)))

    cells = vtk.vtkCellArray()
    conn = np.ascontiguousarray(tris, dtype=np.int64)
    ids = np.column_stack([np.full(len(conn), 3, dtype=np.int64), conn]).ravel()
    cells.SetCells(len(conn), numpy_to_vtk(ids, array_type=vtk.VTK_ID_TYPE))

    poly = vtk.vtkPolyData()
    poly.SetPoints(pts)
    poly.SetPolys(cells)
    tid = numpy_to_vtk(np.ascontiguousarray(tissue_id, dtype=np.int32))
    tid.SetName("tissue_id")
    poly.GetCellData().AddArray(tid)

    # FreeCAD's Fem::FemPostPipeline ViewObject.Field enum only lists POINT
    # data, never cell data -- confirmed live 2026-09-11: loading this file
    # with only the cell array above gets a correct pipeline (right bounds,
    # right tissue_id values via pipe.Data) but Field's enum comes back
    # ['None'] forever, so there is no way to colour by it from the GUI.
    # vtkCellDataToPointData averages each point's neighbouring cells, which
    # for a discrete tissue id blends across boundaries (e.g. grey=9 next to
    # white=10 shows ~9.5 right at the seam) -- fine, even nice, for a visual
    # QA tool; not fine if this array were ever used for anything exact.
    c2p = vtk.vtkCellDataToPointData()
    c2p.SetInputData(poly)
    c2p.Update()
    tid_pt = c2p.GetOutput().GetPointData().GetArray("tissue_id")
    tid_pt.SetName("tissue_id_pt")
    poly.GetPointData().AddArray(tid_pt)

    w = vtk.vtkXMLPolyDataWriter()
    w.SetFileName(path)
    w.SetInputData(poly)
    if w.Write() != 1:
        raise IOError("VTK failed writing " + path)


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: mesh_preview.py <scratch_dir> <params.json>")
    scratch_dir, params_path = sys.argv[1], sys.argv[2]
    if os.path.abspath(scratch_dir) == os.path.abspath(C.OUT):
        raise SystemExit("refusing to use fem/out/ itself as a preview scratch "
                          "dir -- that would overwrite the field of record")

    with open(params_path) as fh:
        raw = json.load(fh)
    params = dict(H_MIN=raw["H_MIN"], H_MAX=raw["H_MAX"],
                  FIELDS=[tuple(x) for x in raw["FIELDS"]],
                  GRID=raw["GRID"], MARGIN=raw["MARGIN"])
    lead_dir = raw.get("LEAD_DIR") or None
    contact_stl, insulator_stl = C.lead_stls(lead_dir)
    if lead_dir is not None and not C.is_single_lead_dir(lead_dir):
        raise SystemExit(
            "LEAD_DIR %r is not a single-lead flat-layout export (missing "
            "its insulator STL, or has zero 'SCS Lead Electrode *.stl' "
            "contact files) -- a multi-lead export directory (per-lead "
            "subdirectories + leads.txt) is out of scope here" % lead_dir)
    # The classification order for THIS lead's own contact count (however
    # many config.lead_stls() found/was told to use above) -- NOT
    # config.ORDER, which is fixed to the frozen 8-contact lead. See
    # config.order_for()'s docstring.
    order = C.order_for(contact_stl)

    def progress(pct):
        print("PROGRESS %d" % pct, flush=True)

    signature = artifacts.mesh_signature(params, contact_stl, insulator_stl)
    progress(0)
    t0 = time.time()
    nodes, tets, lo, hi = build_mesh.build(params=params, out_dir=scratch_dir,
                                           contact_stl=contact_stl,
                                           insulator_stl=insulator_stl,
                                           progress_cb=progress)

    lab = classify_tets(nodes, tets, contact_stl=contact_stl,
                         insulator_stl=insulator_stl, order=order)
    progress(85)
    print("Checking dura barrier (partitioned face sort)", flush=True)

    keep = lab >= 0
    tets_k, lab_k = tets[keep], lab[keep]

    leak = dura_leak_report(nodes, tets_k, lab_k, order=order)
    progress(92)

    print("Extracting tissue interfaces (partitioned face sort)", flush=True)
    tri, tri_lab = tissue_surfaces(nodes, tets_k, lab_k)
    preview_path = os.path.join(scratch_dir, "preview.vtp")
    # Drop interior vertices from the display artifact.
    used, inverse = np.unique(tri, return_inverse=True)
    write_preview_vtu(preview_path + '.tmp.vtp', nodes[used], inverse.reshape(-1, 3), tri_lab)
    os.replace(preview_path + '.tmp.vtp', preview_path)
    progress(100)

    counts = {}
    for i, name in enumerate(order):
        n = int((lab_k == i).sum())
        if n:
            counts[name] = n

    report = dict(
        signature=signature, model_signature=artifacts.model_signature(),
        nodes=len(nodes), tets_total=len(tets), tets_classified=len(tets_k),
        tets_background_dropped=len(tets) - len(tets_k),
        tissue_tet_counts=counts,
        dura=leak,
        preview_vtp=preview_path,
        elapsed_s=round(time.time() - t0, 1),
        params=params,
        lead_dir=lead_dir,
        mesh_npz=os.path.join(scratch_dir, "mesh.npz"),
        # The tissue_id integers embedded in preview.vtp's cell array are
        # positions into THIS list (not config.ORDER) -- SCS_Mesh_Generator.
        # FCMacro's tissue_facet_colors() needs it to colour an N-contact
        # lead's preview correctly, N != 8 included.
        order=order,
        n_contacts=len(contact_stl),
    )
    if signature != artifacts.mesh_signature(params, contact_stl, insulator_stl):
        raise RuntimeError("Model inputs changed during meshing; result not published")
    artifacts.publish(preview_path, signature)
    artifacts.atomic_json(os.path.join(scratch_dir, 'mesh_report.json'), report)
    print("MESH_PREVIEW_JSON " + json.dumps(report))


if __name__ == "__main__":
    main()
