"""Voxelize RADO's spinal cord (white + grey matter, plus the CSF around it)
into NIfTI images that the Spinal Cord Toolbox can register PAM50 against.

    src/atlas/.venv/bin/python src/atlas/voxelize_cord.py [--voxel 0.25]

Writes fem/out/atlas/vox<voxel>/ (git-ignored):
    rado_cord_seg.nii.gz   1 inside white OR grey matter
    rado_gm_seg.nii.gz     1 inside grey matter
    rado_wm_seg.nii.gz     1 inside white matter and not grey
    rado_synth_t2.nii.gz   T2-like contrast: CSF bright, cord dark, GM > WM
    voxelize_report.json   frame, affine, volumes, per-slice shape, laterality

THE FRAME. RADO's STL coordinates are mm in an LPS frame:
    +x = patient LEFT, +y = POSTERIOR (dorsal), +z = SUPERIOR (rostral).
  * +z rostral: the disc files run T_disk_9_10 (z~152) down to T_disk_12_L
    (z~66), and the original fem_dorsal_T10 lead sits on the top vertebra.
  * +y posterior: the aorta and vertebral bodies are at low y, the lead and
    dorsal columns at high y (and see the scs-neuron-dorsal-ventral notes).
  * +x left: (a) the descending thoracic aorta sits LEFT of the vertebral
    midline, and Thoracic_aorta-1.STL's centroid is ~10 mm toward +x of the
    cord; (b) independently, a right-handed CAD frame with +y posterior and
    +z superior forces +x = left. Both are re-checked on every run (the aorta
    test is asserted) because a mirrored atlas would swap the near and far
    dorsal column for a lead 1 mm off the midline, and nothing else would
    visibly fail.
NIfTI world coordinates are RAS, so the affine maps voxel (i, j, k) along
RADO's (x, y, z) to RAS = (-x, -y, z). That is a 180 degree rotation about z,
not a reflection: the image is stored LPS-ordered and every tool that honours
the affine (SCT, FSLeyes, ITK) shows it the right way round. To get back to
RADO mm from RAS world mm: (x, y, z) = (-R, -A, S).

Containment uses fem/scripts/inside.py's InsideTester, one tester per STL
body (assign_and_solve.build_testers explains why bodies must not be merged).
Grey wins over white, matching the FEM classification order. Nothing under
fem/scripts/ is modified; it is only imported.
"""
import argparse
import json
import os
import sys
import time

import nibabel as nib
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "fem", "scripts"))

import config as C  # noqa: E402
from inside import InsideTester  # noqa: E402
from stlio import read_stl  # noqa: E402

AORTA_STL = os.path.join(C.STL, "T8-10 - Thoracic_aorta-1.STL")

# Synthetic T2-weighted intensities. Only the ordering matters to a
# registration cost (CSF >> GM > WM > outside-canal), and PAM50_t2 follows it.
SYNTH = {"outside": 0.0, "csf": 1.0, "gm": 0.45, "wm": 0.25}


def mesh_volume(verts, tris):
    """Enclosed volume of a closed triangle mesh (divergence theorem), mm^3."""
    a, b, c = verts[tris[:, 0]], verts[tris[:, 1]], verts[tris[:, 2]]
    return abs(np.einsum("ij,ij->i", a, np.cross(b, c)).sum()) / 6.0


def inside_any(testers, pts):
    m = np.zeros(len(pts), dtype=bool)
    for t in testers:
        m |= t(pts)
    return m


def load(paths):
    meshes = [read_stl(p) for p in paths]
    return meshes, [InsideTester(v, t) for v, t in meshes]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--voxel", type=float, default=0.25, help="isotropic voxel size, mm")
    ap.add_argument("--margin", type=float, default=1.0,
                    help="padding around the CSF sac bounding box, mm")
    ap.add_argument("--out", default=None, help="output directory")
    args = ap.parse_args()
    h = args.voxel
    out = args.out or os.path.join(C.OUT, "atlas", "vox%g" % h)
    os.makedirs(out, exist_ok=True)
    t0 = time.time()

    white_m, white_t = load(C.TISSUE_BODIES["white"])
    grey_m, grey_t = load(C.TISSUE_BODIES["grey"])
    # The canal's CSF is the main sac only; the other csf bodies are the root
    # sleeves out in the foramina, which are outside this grid anyway.
    csf_m, csf_t = load([C.TISSUE_STL["csf"]])

    # Grid: the CSF sac's bounding box plus a margin, voxel centres on a
    # regular lattice in RADO mm.
    lo = csf_m[0][0].min(0) - args.margin
    hi = csf_m[0][0].max(0) + args.margin
    shape = np.ceil((hi - lo) / h).astype(int)
    axes = [lo[d] + h * (0.5 + np.arange(shape[d])) for d in range(3)]
    X, Y, Z = np.meshgrid(*axes, indexing="ij")
    pts = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    del X, Y, Z
    print("grid %s = %.1f M voxels at %g mm" % (tuple(shape), pts.shape[0] / 1e6, h), flush=True)

    in_grey = inside_any(grey_t, pts)
    in_white = inside_any(white_t, pts) & ~in_grey
    in_csf = inside_any(csf_t, pts) & ~in_grey & ~in_white
    print("containment done in %.0f s" % (time.time() - t0), flush=True)

    gm = in_grey.reshape(shape).astype(np.uint8)
    wm = in_white.reshape(shape).astype(np.uint8)
    cord = gm | wm
    synth = np.full(shape, SYNTH["outside"], dtype=np.float32)
    synth[in_csf.reshape(shape)] = SYNTH["csf"]
    synth[wm.astype(bool)] = SYNTH["wm"]
    synth[gm.astype(bool)] = SYNTH["gm"]

    # voxel (i,j,k) -> RAS world: R = -x, A = -y, S = z.
    affine = np.array([[-h, 0, 0, -axes[0][0]],
                       [0, -h, 0, -axes[1][0]],
                       [0, 0, h, axes[2][0]],
                       [0, 0, 0, 1.0]])

    def save(arr, name):
        img = nib.Nifti1Image(arr, affine)
        img.set_qform(affine, code=1)
        img.set_sform(affine, code=1)
        img.header.set_xyzt_units("mm")
        nib.save(img, os.path.join(out, name))

    save(cord, "rado_cord_seg.nii.gz")
    save(gm, "rado_gm_seg.nii.gz")
    save(wm, "rado_wm_seg.nii.gz")
    save(synth, "rado_synth_t2.nii.gz")

    # --- checks -----------------------------------------------------------
    vv = h ** 3
    vol_white_stl = sum(mesh_volume(*m) for m in white_m)
    vol_grey_stl = sum(mesh_volume(*m) for m in grey_m)
    overlap = int((inside_any(white_t, pts) & in_grey).sum())
    vols = {
        "grey_voxel_mm3": float(gm.sum() * vv),
        "grey_stl_mm3": vol_grey_stl,
        "white_voxel_mm3": float(wm.sum() * vv),
        "white_stl_mm3": vol_white_stl,
        "voxels_inside_both_white_and_grey_shells": overlap,
    }

    # Laterality: the aorta must land on the patient's LEFT, i.e. at larger
    # RADO x (smaller RAS R) than the cord.
    av, at = read_stl(AORTA_STL)
    aorta_x = float(av.mean(0)[0])
    cord_x = float(pts[cord.ravel().astype(bool), 0].mean())
    ras_aorta = affine @ np.r_[(np.array([aorta_x, 0, 0]) - [axes[0][0], axes[1][0], axes[2][0]]) / h, 1]
    ras_cord = affine @ np.r_[(np.array([cord_x, 0, 0]) - [axes[0][0], axes[1][0], axes[2][0]]) / h, 1]
    lat = {
        "aorta_centroid_x_mm": aorta_x,
        "cord_centroid_x_mm": cord_x,
        "aorta_minus_cord_x_mm": aorta_x - cord_x,
        "aorta_R_minus_cord_R_in_RAS_mm": float(ras_aorta[0] - ras_cord[0]),
        "verdict": "aorta is patient-left of the cord, as expected" if aorta_x > cord_x
                   else "MIRRORED: aorta is patient-right; do not trust this frame",
    }
    assert aorta_x > cord_x + 3.0, lat["verdict"]
    assert nib.aff2axcodes(affine) == ("L", "P", "S")

    # Per-slice cord shape, for matching against PAM50 levels later.
    slices = []
    for k in range(0, shape[2], max(1, int(round(1.0 / h)))):
        c = cord[:, :, k].astype(bool)
        if not c.any():
            continue
        g = gm[:, :, k].astype(bool)
        ii, jj = np.nonzero(c)
        slices.append({
            "z_mm": round(float(axes[2][k]), 3),
            "cord_area_mm2": round(float(c.sum() * h * h), 3),
            "gm_area_mm2": round(float(g.sum() * h * h), 3),
            "lr_width_mm": round(float((ii.max() - ii.min() + 1) * h), 3),
            "ap_depth_mm": round(float((jj.max() - jj.min() + 1) * h), 3),
            "centroid_x_mm": round(float(axes[0][ii].mean()), 3),
            "centroid_y_mm": round(float(axes[1][jj].mean()), 3),
        })

    report = {
        "voxel_mm": h,
        "shape": shape.tolist(),
        "rado_origin_of_voxel_000_mm": [float(a[0]) for a in axes],
        "frame": "RADO mm is LPS (+x left, +y posterior, +z superior); "
                 "NIfTI world RAS = (-x, -y, z)",
        "affine_vox_to_ras": affine.tolist(),
        "axcodes": list(nib.aff2axcodes(affine)),
        "sources": {"white": C.TISSUE_BODIES["white"], "grey": C.TISSUE_BODIES["grey"],
                    "csf": [C.TISSUE_STL["csf"]], "laterality_landmark": AORTA_STL},
        "volumes": vols,
        "laterality": lat,
        "per_slice_every_1mm": slices,
        "elapsed_s": round(time.time() - t0, 1),
    }
    with open(os.path.join(out, "voxelize_report.json"), "w") as f:
        json.dump(report, f, indent=1)

    print("grey  %.1f mm3 voxel vs %.1f mm3 STL (%+.2f%%)" % (
        vols["grey_voxel_mm3"], vol_grey_stl, 100 * (vols["grey_voxel_mm3"] / vol_grey_stl - 1)))
    print("white %.1f mm3 voxel vs %.1f mm3 STL (white shell also counted: %+.2f%%)" % (
        vols["white_voxel_mm3"], vol_white_stl,
        100 * ((vols["white_voxel_mm3"] + overlap * vv) / vol_white_stl - 1)))
    print("voxels inside both shells: %d" % overlap)
    print("laterality: %s (aorta %+.1f mm in x from cord)" % (lat["verdict"], aorta_x - cord_x))
    print("wrote %s in %.0f s" % (out, time.time() - t0))


if __name__ == "__main__":
    main()
