"""Visual check of voxelize_cord.py output: axial slices at each poster lead's
driven contacts, drawn in patient coordinates read back from the NIfTI affine
(not from RADO's axes), with that run's contacts overlaid.

    src/atlas/.venv/bin/python src/atlas/plot_voxel_check.py [--voxel 0.25]

If the frame were mirrored, the contacts would appear on the wrong side of the
midline relative to their RADO x, and the 'patient L' label would sit on the
aorta-free side. Writes fem/out/atlas/vox<voxel>/voxel_check.png.
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

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "fem", "scripts"))

import config as C  # noqa: E402
from stlio import read_stl  # noqa: E402

# Poster runs (historical provenance; geometry only is read here) and the
# z of the cross-sections used on the poster.
MIDLINE_X = 56.35  # cord centroid x, RADO mm (voxelize_report.json)
RUNS = {"dorsal": ("63c73a9ab23d01d11669c7b4", 78.5),
        "ventral": ("1b4c9eec2c617fe8e4404f4a", 136.8)}


def contacts(run):
    lead = os.path.join(C.LEAD_RUNS, run, "lead")
    out = {}
    for n in range(1, 17):
        p = os.path.join(lead, "SCS Lead Electrode %d.stl" % n)
        if os.path.exists(p):
            out[n] = read_stl(p)[0].mean(0)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--voxel", type=float, default=0.25)
    args = ap.parse_args()
    d = os.path.join(C.OUT, "atlas", "vox%g" % args.voxel)
    img = nib.load(os.path.join(d, "rado_synth_t2.nii.gz"))
    vol = np.asarray(img.dataobj)
    A = img.affine
    inv = np.linalg.inv(A)
    with open(os.path.join(d, "disc_labels.json")) as f:
        discs = json.load(f)["labels"]

    fig, axs = plt.subplots(1, 3, figsize=(15, 5.6),
                            gridspec_kw={"width_ratios": [1, 1, 0.9]})
    for ax in axs:
        ax.set_facecolor("0.35")  # grey outside the voxel grid, not white
    for ax, (name, (run, z)) in zip(axs[:2], RUNS.items()):
        k = int(round((inv @ [0, 0, z, 1])[2]))
        sl = vol[:, :, k]
        # world RAS of each voxel centre in this slice
        ni, nj = sl.shape
        corner = [A @ [i, j, k, 1] for i, j in ((-.5, -.5), (ni - .5, nj - .5))]
        R0, A0 = corner[0][0], corner[0][1]
        R1, A1 = corner[1][0], corner[1][1]
        # Show in RADIOLOGICAL convention: patient RIGHT on the viewer's LEFT,
        # anterior at the bottom. imshow wants [row=A, col=R].
        ax.imshow(sl.T, origin="lower", cmap="gray", extent=(R0, R1, A0, A1),
                  interpolation="nearest")
        ax.set_xlim(max(R0, R1) + 6, min(R0, R1) - 6)  # +R on the left
        ax.set_ylim(min(A0, A1) - 6, max(A0, A1) + 8)
        ax.invert_yaxis()                               # anterior at bottom
        for n, c in contacts(run).items():
            if abs(c[2] - z) > 6:
                continue
            ras = np.array([-c[0], -c[1], c[2]])
            ax.plot(ras[0], ras[1], "o", ms=9, mfc="none",
                    mec="tab:red" if n <= 8 else "tab:cyan", mew=2)
            ax.annotate("C%d" % n, (ras[0], ras[1]), xytext=(4, 4),
                        textcoords="offset points", color="w", fontsize=8)
        ax.axvline(-MIDLINE_X, color="y", lw=0.6, ls="--")
        ax.set_title("%s run, z = %.1f mm (RADO)" % (name, z))
        ax.set_xlabel("RAS R (mm)   <- patient RIGHT | patient LEFT ->")
        ax.set_ylabel("RAS A (mm)   anterior at bottom")
        ax.text(0.02, 0.97, "patient R", transform=ax.transAxes, color="w", va="top")
        ax.text(0.98, 0.97, "patient L", transform=ax.transAxes, color="w", va="top",
                ha="right")

    # sagittal through the cord midline, with disc labels
    ax = axs[2]
    i = int(round((inv @ [-MIDLINE_X, 0, 0, 1])[0]))
    sl = vol[i, :, :]
    nj, nk = sl.shape
    c0, c1 = A @ [i, -.5, -.5, 1], A @ [i, nj - .5, nk - .5, 1]
    ax.imshow(sl.T, origin="lower", cmap="gray", extent=(c0[1], c1[1], c0[2], c1[2]),
              aspect="equal", interpolation="nearest")
    ax.set_xlim(min(c0[1], c1[1]), max(c0[1], c1[1]))  # anterior to the right
    for e in discs:
        x, y, zz = e["label_rado_mm"]
        ax.plot(-y, zz, "r+", ms=10, mew=2)
        ax.annotate("%s (%d)" % (e["vertebra_below_by_disc_name"], e["label_value"]),
                    (-y, zz), xytext=(6, -3), textcoords="offset points", color="r", fontsize=8)
    for name, (run, z) in RUNS.items():
        ax.axhline(z, color="tab:orange", lw=0.8, ls=":")
        ax.text(ax.get_xlim()[0], z + 1, name + " slice", color="tab:orange", fontsize=8)
    ax.set_title("midsagittal; red + = disc labels (vertebra below)")
    ax.set_xlabel("RAS A (mm)   posterior <- | -> anterior")
    ax.set_ylabel("S (mm)")

    fig.suptitle("RADO cord voxelized (synthetic T2: CSF bright, GM > WM); "
                 "patient sides from the NIfTI affine; red = lead 1 contacts, cyan = lead 2",
                 fontsize=10)
    fig.tight_layout()
    out = os.path.join(d, "voxel_check.png")
    fig.savefig(out, dpi=130)
    print(out)


if __name__ == "__main__":
    main()
