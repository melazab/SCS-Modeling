"""Which PAM50 levels have RADO's cord cross-section? A shape screen to inform
the vertebral-level-versus-shape choice for template registration.

    src/atlas/.venv/bin/python src/atlas/pam50_level_match.py

RADO's cord is a single cross-section swept along the canal (voxelize_cord.py
per-slice report: 8.25 x 5 mm, 32 mm2 cord, 9.6 mm2 grey at every z). So "which
PAM50 slice looks like RADO" is one comparison per PAM50 slice, not per z.

For each PAM50 axial slice (PAM50 is straightened, so its slices are
perpendicular to the cord):
  1. Normalise both cords to their own centroid and left-right / antero-
     posterior half-extents, so overall size and aspect ratio drop out.
  2. Sample both on the same normalised grid; grey = PAM50_gm >= 0.5.
  3. Report grey-matter Dice and cord-outline Dice in that normalised frame,
     plus raw areas, grey fraction and aspect ratio.
RADO is sampled straight from the STLs (InsideTester) on a plane
PERPENDICULAR to its centreline at --z (centreline from voxelize_cord.py,
which must have been run at 0.25 mm), because RADO's z-slices are ~12 deg
oblique to the cord.

This is a screen, not a registration: the Dice that SCT achieves after its own
non-linear slice-wise warp will be higher everywhere. What matters here is the
ranking along the cord. Writes fem/out/atlas/pam50_level_match.{json,png}.
"""
import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import nibabel as nib  # noqa: E402
import numpy as np  # noqa: E402
from scipy.ndimage import map_coordinates  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "fem", "scripts"))

import config as C  # noqa: E402
from inside import InsideTester  # noqa: E402
from stlio import read_stl  # noqa: E402

PAM50 = os.path.expanduser("~/sct_7.3/data/PAM50/template")
VERT = {**{i: "C%d" % i for i in range(1, 8)}, **{7 + i: "T%d" % i for i in range(1, 13)},
        **{19 + i: "L%d" % i for i in range(1, 6)}}
SPINAL = {**{i: "C%d" % i for i in range(1, 9)}, **{8 + i: "T%d" % i for i in range(1, 13)},
          **{20 + i: "L%d" % i for i in range(1, 6)}, **{25 + i: "S%d" % i for i in range(1, 6)}}
N = 121            # normalised grid points per axis
SPAN = 1.15        # grid covers +/- SPAN half-extents


def dice(a, b):
    s = a.sum() + b.sum()
    return float(2 * (a & b).sum() / s) if s else np.nan


def rado_section(z0, h=0.05):
    """Cord and grey masks of RADO on the plane through the centreline at z0,
    perpendicular to it. Returns masks on a (left, posterior) mm grid."""
    white = [InsideTester(*read_stl(p)) for p in C.TISSUE_BODIES["white"]]
    grey = [InsideTester(*read_stl(p)) for p in C.TISSUE_BODIES["grey"]]
    # Centreline from voxelize_cord.py's per-slice cord centroids. (Not from
    # STL vertices: the swept STLs use long triangles, so mid-cord z-bands can
    # hold no vertices at all.)
    with open(os.path.join(C.OUT, "atlas", "vox0.25", "voxelize_report.json")) as f:
        sl = json.load(f)["per_slice_every_1mm"]
    zs = np.array([e["z_mm"] for e in sl])
    cx = np.array([e["centroid_x_mm"] for e in sl])
    cy = np.array([e["centroid_y_mm"] for e in sl])

    def centre(z):
        return np.array([np.interp(z, zs, cx), np.interp(z, zs, cy), z])
    c0 = centre(z0)
    t = centre(z0 + 3) - centre(z0 - 3)
    t /= np.linalg.norm(t)
    ex = np.array([1.0, 0, 0])                      # patient left, unchanged
    ey = np.cross(t, ex)                            # in-plane, ~ +y (posterior)
    ey /= np.linalg.norm(ey)
    if ey[1] < 0:
        ey = -ey
    u = np.arange(-6, 6 + h / 2, h)
    U, V = np.meshgrid(u, u, indexing="ij")
    pts = c0 + U.ravel()[:, None] * ex + V.ravel()[:, None] * ey
    g = np.zeros(len(pts), bool)
    w = np.zeros(len(pts), bool)
    for tt in grey:
        g |= tt(pts)
    for tt in white:
        w |= tt(pts)
    tilt = np.degrees(np.arccos(abs(t[2])))
    return u, (g | w).reshape(U.shape), g.reshape(U.shape), float(tilt), c0


def normalise_grid(cord, coords_l, coords_p):
    """Centroid and half-extents of a cord mask given per-axis mm coords."""
    ii, jj = np.nonzero(cord)
    cl, cp = coords_l[ii].mean(), coords_p[jj].mean()
    hl = (coords_l[ii].max() - coords_l[ii].min()) / 2
    hp = (coords_p[jj].max() - coords_p[jj].min()) / 2
    return cl, cp, hl, hp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--z", type=float, default=110.0, help="RADO z of the reference section, mm")
    args = ap.parse_args()
    out = os.path.join(C.OUT, "atlas")
    os.makedirs(out, exist_ok=True)

    # --- RADO reference section --------------------------------------------
    u, rc, rg, tilt, c0 = rado_section(args.z)
    h = u[1] - u[0]
    cl, cp, hl, hp = normalise_grid(rc, u, u)
    s = np.linspace(-SPAN, SPAN, N)
    SL, SP = np.meshgrid(s, s, indexing="ij")
    # nearest sample of the fine RADO grid at the normalised points
    il = np.clip(np.round((cl + SL * hl - u[0]) / h).astype(int), 0, len(u) - 1)
    ip = np.clip(np.round((cp + SP * hp - u[0]) / h).astype(int), 0, len(u) - 1)
    R_cord, R_gm = rc[il, ip], rg[il, ip]
    rado = {"z_mm": args.z, "centreline_point_mm": c0.round(3).tolist(),
            "tilt_from_z_deg": round(tilt, 2),
            "cord_area_mm2": float(rc.sum() * h * h), "gm_area_mm2": float(rg.sum() * h * h),
            "lr_width_mm": 2 * hl, "ap_depth_mm": 2 * hp}
    rado["gm_fraction"] = rado["gm_area_mm2"] / rado["cord_area_mm2"]
    rado["aspect_lr_over_ap"] = hl / hp

    # --- PAM50, slice by slice ---------------------------------------------
    cord_img = nib.load(os.path.join(PAM50, "PAM50_cord.nii.gz"))
    assert nib.aff2axcodes(cord_img.affine) == ("L", "A", "S")
    cord = np.asarray(cord_img.dataobj) > 0
    gm = np.asarray(nib.load(os.path.join(PAM50, "PAM50_gm.nii.gz")).dataobj)
    lev = np.asarray(nib.load(os.path.join(PAM50, "PAM50_levels.nii.gz")).dataobj)
    spl = np.asarray(nib.load(os.path.join(PAM50, "PAM50_spinal_levels.nii.gz")).dataobj)
    vx = float(cord_img.header.get_zooms()[0])
    nx, ny, nz = cord.shape
    # LAS voxel axes: i increases toward patient LEFT, j toward ANTERIOR.
    # Use (left, posterior) mm coordinates to match RADO's frame.
    coords_l = np.arange(nx) * vx
    coords_p = -np.arange(ny) * vx
    rows = []
    for k in range(nz):
        c = cord[:, :, k]
        if c.sum() < 20:
            continue
        L0, P0, HL, HP = normalise_grid(c, coords_l, coords_p)
        fi = (L0 + SL * HL) / vx
        fj = -(P0 + SP * HP) / vx
        Pc = map_coordinates(c.astype(np.float32), [fi, fj], order=1) >= 0.5
        Pg = map_coordinates(gm[:, :, k], [fi, fj], order=1) >= 0.5
        vl = lev[:, :, k][c]
        sl = spl[:, :, k][c]
        vl = int(np.bincount(vl).argmax()) if vl.any() else 0
        sl = int(np.bincount(sl).argmax()) if sl.any() else 0
        rows.append({"k": k, "vert": VERT.get(vl, ""), "spinal": SPINAL.get(sl, ""),
                     "cord_area_mm2": float(c.sum() * vx * vx),
                     "gm_area_mm2": float(gm[:, :, k][c].sum() * vx * vx),
                     "aspect_lr_over_ap": HL / HP,
                     "gm_dice_normalised": dice(Pg, R_gm),
                     "cord_dice_normalised": dice(Pc, R_cord)})
        rows[-1]["gm_fraction"] = rows[-1]["gm_area_mm2"] / rows[-1]["cord_area_mm2"]

    # summary per vertebral level and per spinal level
    def summarise(key):
        out_ = {}
        for r in rows:
            if r[key]:
                out_.setdefault(r[key], []).append(r)
        return {lv: {"n_slices": len(v),
                     "gm_dice_median": float(np.median([r["gm_dice_normalised"] for r in v])),
                     "gm_fraction_median": float(np.median([r["gm_fraction"] for r in v])),
                     "cord_area_median_mm2": float(np.median([r["cord_area_mm2"] for r in v])),
                     "aspect_median": float(np.median([r["aspect_lr_over_ap"] for r in v]))}
                for lv, v in out_.items()}

    rep = {"rado_reference": rado, "by_vertebral_level": summarise("vert"),
           "by_spinal_level": summarise("spinal"), "slices": rows,
           "method": __doc__}
    with open(os.path.join(out, "pam50_level_match.json"), "w") as f:
        json.dump(rep, f, indent=1)

    # --- figure --------------------------------------------------------------
    ks = np.array([r["k"] for r in rows])
    fig, axs = plt.subplots(4, 1, figsize=(11, 10), sharex=True)
    for ax, key, lab, ref in (
            (axs[0], "gm_dice_normalised", "grey Dice vs RADO\n(normalised)", None),
            (axs[1], "gm_fraction", "grey / cord area", rado["gm_fraction"]),
            (axs[2], "cord_area_mm2", "cord area, mm2", rado["cord_area_mm2"]),
            (axs[3], "aspect_lr_over_ap", "LR / AP", rado["aspect_lr_over_ap"])):
        ax.plot(ks, [r[key] for r in rows], lw=1)
        if ref is not None:
            ax.axhline(ref, color="tab:red", ls="--", lw=1, label="RADO")
            ax.legend(loc="upper left", fontsize=8)
        ax.set_ylabel(lab)
        ax.grid(alpha=0.3)
    # level ticks
    for key, y, col in (("vert", 1.02, "k"), ("spinal", 1.10, "tab:blue")):
        seen = {}
        for r in rows:
            seen.setdefault(r[key], []).append(r["k"])
        for lv, kk in seen.items():
            if lv:
                axs[0].text(np.mean(kk), y, lv, transform=axs[0].get_xaxis_transform(),
                            ha="center", fontsize=6, color=col)
                axs[0].axvline(min(kk), color=col, lw=0.3, alpha=0.4)
    axs[0].set_title("PAM50 vs RADO cross-section; black = vertebral level, blue = spinal "
                     "(cord) segment\n", fontsize=10)
    axs[-1].set_xlabel("PAM50 slice (0.5 mm, caudal <- -> rostral)")
    axs[-1].set_xlim(ks.min(), ks.max())
    axs[-1].invert_xaxis()
    fig.tight_layout()
    fig.savefig(os.path.join(out, "pam50_level_match.png"), dpi=130)

    print("RADO section at z=%.0f (tilt %.1f deg): cord %.1f mm2, grey %.1f mm2 (%.2f), "
          "%.2f x %.2f mm" % (args.z, tilt, rado["cord_area_mm2"], rado["gm_area_mm2"],
                              rado["gm_fraction"], rado["lr_width_mm"], rado["ap_depth_mm"]))
    print("\nby vertebral level:   grey Dice   grey frac   cord mm2   LR/AP")
    for lv, v in rep["by_vertebral_level"].items():
        print("  %-4s %14.3f %11.3f %10.1f %7.2f" % (lv, v["gm_dice_median"], v["gm_fraction_median"],
                                                     v["cord_area_median_mm2"], v["aspect_median"]))
    print("\nby spinal segment:")
    for lv, v in rep["by_spinal_level"].items():
        print("  %-4s %14.3f %11.3f %10.1f %7.2f" % (lv, v["gm_dice_median"], v["gm_fraction_median"],
                                                     v["cord_area_median_mm2"], v["aspect_median"]))


if __name__ == "__main__":
    main()
