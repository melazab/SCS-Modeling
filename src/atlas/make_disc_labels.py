"""Disc labels for SCT template registration, on voxelize_cord.py's grid.

    src/atlas/.venv/bin/python src/atlas/make_disc_labels.py [--voxel 0.25]
                                                             [--naming discs|paper]

Writes fem/out/atlas/vox<voxel>/rado_disc_labels.nii.gz: one voxel per
intervertebral disc, at the orthogonal projection of the disc onto the cord
(the centre of the cord on the slice through the disc's posterior tip), which
is one of the two placements sct_register_to_template -ldisc accepts.

SCT's value convention: a disc is labelled with the index of the vertebra BELOW
it, counting C1 = 1 (C2/C3 -> 3, T1/T2 -> 9, so T9/T10 -> 17, T12/L1 -> 20).

LEVEL NAMES ARE NOT SETTLED. RADO names its discs T_disk_9_10 ... T_disk_12_L,
which makes the three vertebrae T10-T12; the paper says T9-T11; the STL
filenames say "T8-10". src/freecad/README.md "Vertebral levels" records the
same three-way disagreement and the Lead Designer follows the disc names, so
the default here does too (--naming discs). --naming paper shifts every label
up one level (T9-T11). The z positions are geometry and do not change.
"""
import argparse
import json
import os
import re
import sys

import nibabel as nib
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "fem", "scripts"))

import config as C  # noqa: E402
from stlio import read_stl  # noqa: E402

# C1..C7 = 1..7, T1..T12 = 8..19, L1..L5 = 20..24
VERT_INDEX = {**{"C%d" % i: i for i in range(1, 8)},
              **{"T%d" % i: 7 + i for i in range(1, 13)},
              **{"L%d" % i: 19 + i for i in range(1, 6)}}


def disc_below_vertebra(name):
    """'T_disk_9_10x-1' -> 'T10'; 'T_disk_12_Lx-1' -> 'L1'."""
    m = re.search(r"T_disk_(\d+)_(\d+|L)", name)
    lo = m.group(2)
    return "L1" if lo == "L" else "T%s" % lo


def posterior_tip(verts, midline_x, halfwidth=1.5, depth=0.5):
    """Mean of the disc vertices within `halfwidth` mm of the midline and
    within `depth` mm of the disc's most posterior (max-y) extent there."""
    mid = verts[np.abs(verts[:, 0] - midline_x) < halfwidth]
    return mid[mid[:, 1] > mid[:, 1].max() - depth].mean(0)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--voxel", type=float, default=0.25)
    ap.add_argument("--naming", choices=["discs", "paper"], default="discs")
    ap.add_argument("--dir", default=None, help="voxelize_cord.py output directory")
    args = ap.parse_args()
    d = args.dir or os.path.join(C.OUT, "atlas", "vox%g" % args.voxel)

    cord_img = nib.load(os.path.join(d, "rado_cord_seg.nii.gz"))
    cord = np.asarray(cord_img.dataobj).astype(bool)
    with open(os.path.join(d, "voxelize_report.json")) as f:
        rep = json.load(f)
    h = rep["voxel_mm"]
    o = np.array(rep["rado_origin_of_voxel_000_mm"])
    midline_x = o[0] + h * np.nonzero(cord)[0].mean()

    labels = np.zeros(cord.shape, dtype=np.uint8)
    out = []
    for p in C.TISSUE_BODIES["disc"]:
        name = os.path.basename(p)
        v, _ = read_stl(p)
        tip = posterior_tip(v, midline_x)
        k = int(round((tip[2] - o[2]) / h))
        ii, jj = np.nonzero(cord[:, :, k])
        assert ii.size, "no cord on the slice through %s" % name
        i, j = int(round(ii.mean())), int(round(jj.mean()))
        below = disc_below_vertebra(name)
        value = VERT_INDEX[below] - (1 if args.naming == "paper" else 0)
        labels[i, j, k] = value
        out.append({"stl": name, "vertebra_below_by_disc_name": below,
                    "label_value": value,
                    "posterior_tip_rado_mm": tip.round(3).tolist(),
                    "label_voxel": [i, j, k],
                    "label_rado_mm": (o + h * np.array([i, j, k])).round(3).tolist()})

    img = nib.Nifti1Image(labels, cord_img.affine, cord_img.header)
    nib.save(img, os.path.join(d, "rado_disc_labels.nii.gz"))
    out.sort(key=lambda e: -e["label_rado_mm"][2])
    with open(os.path.join(d, "disc_labels.json"), "w") as f:
        json.dump({"naming": args.naming, "sct_convention": "value = index of vertebra "
                   "below the disc, C1 = 1", "labels": out}, f, indent=1)
    for e in out:
        print("%-20s value %2d  z %.2f mm" % (e["stl"][:-4], e["label_value"], e["label_rado_mm"][2]))


if __name__ == "__main__":
    main()
