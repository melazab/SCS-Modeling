"""Register the PAM50 white-matter atlas onto RADO's cord cross-section.

    src/atlas/.venv/bin/python src/atlas/register_tracts.py [--level T10] [--s 70]

RADO's cord is one cross-section swept along the canal (README "RADO's cord is
one cross-section, swept"), and the extended T6-T12 model copies it. So one 2D
registration covers the whole cord: this script registers ONE PAM50 level onto
that section and the result can be swept along the centreline. Which PAM50
level is the --level choice (README "Consequence for registration"); the
default T10 is the vertebral level of the leads.

Steps:
  1. RADO section: white/grey STLs (RADO's own) sampled on the plane
     perpendicular to the extended model's centreline (FIBRE_CENTRELINE) at
     arc length --s, in that centreline's frame: u = +x (patient left),
     v = t x u (posterior). The centreline lies in x = x0 and the cord is
     swept along it, so these (u, v) place the section at any arc length.
  2. PAM50 section: the middle slice of the chosen vertebral level (PAM50 is
     straightened, so its slices are perpendicular to the cord), its cord,
     grey and the 36 atlas maps, moved into the same (left, posterior) frame.
  3. Both become a synthetic image (white = 1, grey = 2, outside = 0) on the
     same 0.05 mm grid, centroids aligned. SCT registers PAM50 -> RADO:
     affine, then SyN, MeanSquares, on images blurred by BLUR_MM (no real contrast exists, so only
     the white/grey outlines are matched; tract layout inside is the atlas'
     prior). Slabs are 5 identical slices (ITK smoothing needs >= 4); SCT's 2D
     algorithms only deform in-plane.
  4. The 36 atlas maps are warped with the same field (linear), masked to
     RADO's white (tracts 0-29) and grey (30-35), and the white-matter tract
     probabilities renormalised to sum to 1 inside RADO's white matter (the
     atlas has partial volume at the edge). Pixels where the warped atlas has
     no tract at all take the nearest pixel's values.

Writes fem/out/atlas/tracts_<level>/:
  tract_map.npz      u, v (mm), prob [36, nu, nv], label (argmax; -1 outside)
  tracts.json        tract names, colours, areas, Dice before/after, and each
                     tract's argmax region as polygons in (u, v) mm; read by
                     src/freecad/SCS_TractVisualizer.FCMacro
  tract_map.png      PAM50 source, RADO target, warped labels
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import nibabel as nib  # noqa: E402
import numpy as np  # noqa: E402
from scipy.ndimage import distance_transform_edt, gaussian_filter, map_coordinates  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "fem", "scripts"))
os.environ["SCS_GEOMETRY"] = "rado"
import config as C  # noqa: E402
from inside import InsideTester  # noqa: E402
from stlio import read_stl  # noqa: E402

SCT = os.path.expanduser("~/sct_7.3/bin")
PAM50 = os.path.expanduser("~/sct_7.3/data/PAM50")
VERT = {7 + i: "T%d" % i for i in range(1, 13)}
H = 0.05                       # grid, mm
HALF = 6.0                     # grid covers +/- HALF mm
BLUR_MM = 0.25
# Swept on T10, grey Dice after: affine + SyN 0.775 at blur 0.25 mm, against
# 0.70-0.75 for BSplineSyN or blur 0.1 / 0.4 and 0.40 with no blur.
PARAM = ("step=1,type=im,algo=affine,metric=MeanSquares:"
         "step=2,type=im,algo=syn,metric=MeanSquares,iter=100,shrink=1,smooth=0,gradStep=0.2")
# One colour per tract pair (left and right share it); grey regions in greys.
COLOURS = {"fasciculus gracilis": "#1f5fbf", "fasciculus cuneatus": "#57a8e8",
           "lateral corticospinal tract": "#d62728", "ventral spinocerebellar tract": "#ff9f1c",
           "rubrospinal tract": "#8c564b", "lateral reticulospinal tract": "#c49c94",
           "spinal lemniscus (spinothalamic and spinoreticular tracts)": "#2ca02c",
           "spino-olivary tract": "#98df8a", "ventrolateral reticulospinal tract": "#bcbd22",
           "lateral vestibulospinal tract": "#9467bd", "ventral reticulospinal tract": "#c5b0d5",
           "ventral corticospinal tract": "#e377c2", "tectospinal tract": "#17becf",
           "medial reticulospinal tract": "#9edae5", "medial longitudinal fasciculus": "#7f7f7f",
           "ventral horn": "#555555", "intermediate zone": "#888888", "dorsal horn": "#aaaaaa"}


def atlas_names():
    names = []
    with open(os.path.join(PAM50, "atlas", "info_label.txt")) as f:
        for line in f:
            p = [x.strip() for x in line.split(",")]
            if len(p) == 3 and p[0].isdigit() and int(p[0]) < 36:
                names.append((int(p[0]), p[1], p[2]))
    return names


def frame_at(s):
    with open(C.FIBRE_CENTRELINE) as f:
        cl = json.load(f)
    S, P, T = (np.asarray(cl[k], float) for k in ("s_mm", "points_mm", "tangents"))
    c = np.array([np.interp(s, S, P[:, d]) for d in range(3)])
    t = np.array([np.interp(s, S, T[:, d]) for d in range(3)])
    t /= np.linalg.norm(t)
    u = np.array([1.0, 0, 0]) - t[0]*t
    u /= np.linalg.norm(u)
    v = np.cross(t, u)
    assert v[1] > 0, "v should point posterior"
    return c, u, v, t


def rado_section(s, g):
    c, u, v, _ = frame_at(s)
    U, V = np.meshgrid(g, g, indexing="ij")
    pts = c + U.ravel()[:, None]*u + V.ravel()[:, None]*v
    inside = lambda key: np.any([InsideTester(*read_stl(p))(pts) for p in C.TISSUE_BODIES[key]], axis=0)
    grey = inside("grey").reshape(U.shape)
    white = inside("white").reshape(U.shape) & ~grey
    return white, grey, c


def pam50_section(level, g):
    lev = nib.load(os.path.join(PAM50, "template", "PAM50_levels.nii.gz"))
    assert nib.aff2axcodes(lev.affine) == ("L", "A", "S")
    L = np.asarray(lev.dataobj)
    code = [k for k, n in VERT.items() if n == level][0]
    ks = np.nonzero((L == code).any(axis=(0, 1)))[0]
    k = int(np.median(ks))
    vx = float(lev.header.get_zooms()[0])
    load = lambda p: np.asarray(nib.load(p).dataobj)[:, :, k].astype(np.float32)
    cord = load(os.path.join(PAM50, "template", "PAM50_cord.nii.gz"))
    ii, jj = np.nonzero(cord > 0.5)
    ci, cj = ii.mean(), jj.mean()
    # LAS voxels: i -> left, j -> anterior. Output grid is (left, posterior).
    U, V = np.meshgrid(g, g, indexing="ij")
    fi, fj = ci + U/vx, cj - V/vx

    def sample(img):
        return map_coordinates(img, [fi, fj], order=1)
    maps = [sample(load(os.path.join(PAM50, "atlas", fn))) for _, _, fn in atlas_names()]
    return sample(cord), sample(load(os.path.join(PAM50, "template", "PAM50_gm.nii.gz"))), np.array(maps), k


def to_nii(a, path):
    """2D (u, v) array as a 5-slice slab. Voxel axes (u = left, v = posterior,
    slab); world RAS = (-u, -v, slab)."""
    vol = np.repeat(a[:, :, None], 5, axis=2).astype(np.float32)
    aff = np.diag([-H, -H, H, 1.0])
    aff[:2, 3] = [HALF, HALF]       # u = g[0] = -HALF at voxel 0 -> RAS x = +HALF
    nib.save(nib.Nifti1Image(vol, aff), path)


def from_nii(path):
    return np.asarray(nib.load(path).dataobj)[:, :, 2]


def dice(a, b):
    return float(2*(a & b).sum()/(a.sum() + b.sum()))


def polygons(mask):
    """Outer and hole rings of a boolean mask on the (u, v) grid, in mm."""
    fig = plt.figure()
    g = np.arange(-HALF, HALF + H/2, H)
    cs = plt.contour(g, g, mask.T.astype(float), levels=[0.5])
    rings = [p.tolist() for p in cs.allsegs[0] if len(p) >= 4]
    plt.close(fig)
    return [[[round(x, 4), round(y, 4)] for x, y in r] for r in rings]


def run(cmd):
    print("+", " ".join(cmd))
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", default="T10", help="PAM50 vertebral level used (T1-T12)")
    ap.add_argument("--s", type=float, default=70.0, help="extended-centreline arc length of the RADO section, mm")
    args = ap.parse_args()
    out = os.path.join(C.OUT, "atlas", "tracts_%s" % args.level)
    os.makedirs(out, exist_ok=True)
    g = np.arange(-HALF, HALF + H/2, H)

    rw, rg, c0 = rado_section(args.s, g)
    rcord = rw | rg
    # centre of the RADO cord in the centreline frame
    U, V = np.meshgrid(g, g, indexing="ij")
    offset = np.array([U[rcord].mean(), V[rcord].mean()])
    pc, pg, maps, k = pam50_section(args.level, g)
    # move PAM50 so its cord centroid sits on RADO's
    shift = lambda a: map_coordinates(a, [(U - offset[0] - g[0])/H, (V - offset[1] - g[0])/H], order=1)
    pc, pg = shift(pc), shift(pg)
    maps = np.array([shift(m) for m in maps])
    src = (pc > 0.5).astype(np.float32) + (pg > 0.5)
    dst = rw.astype(np.float32) + 2*rg

    with tempfile.TemporaryDirectory(dir=os.environ.get("SCS_TMP")) as tmp:
        P = lambda n: os.path.join(tmp, n)
        # Blurred label images give the metric a basin wider than the
        # 0.05 mm edge band, so large grey-matter moves can converge.
        blur = float(os.environ.get("SCS_BLUR_MM", BLUR_MM))/H
        to_nii(gaussian_filter(src, blur), P("src.nii.gz")); to_nii(gaussian_filter(dst, blur), P("dst.nii.gz"))
        run([os.path.join(SCT, "sct_register_multimodal"), "-i", P("src.nii.gz"), "-d", P("dst.nii.gz"),
             "-param", os.environ.get("SCS_SCT_PARAM", PARAM), "-ofolder", tmp, "-owarp", P("warp.nii.gz"), "-x", "linear", "-v", os.environ.get("SCT_V", "0")])
        warped = []
        for n, m in enumerate(list(maps) + [pg]):
            to_nii(m, P("m.nii.gz"))
            run([os.path.join(SCT, "sct_apply_transfo"), "-i", P("m.nii.gz"), "-d", P("dst.nii.gz"),
                 "-w", P("warp.nii.gz"), "-x", "linear", "-o", P("m_reg.nii.gz"), "-v", "0"])
            warped.append(from_nii(P("m_reg.nii.gz")))
    wg = warped.pop()
    prob = np.clip(np.array(warped), 0, None)

    names = atlas_names()
    wm_ids, gm_ids = list(range(30)), list(range(30, 36))
    for ids, mask in ((wm_ids, rw), (gm_ids, rg)):
        p = prob[ids]*mask
        tot = p.sum(0)
        empty = mask & (tot < 1e-3)
        if empty.any():                       # nearest pixel that has the atlas
            _, (ni, nj) = distance_transform_edt(~(mask & (tot >= 1e-3)), return_indices=True)
            p[:, empty] = p[:, ni[empty], nj[empty]]
            tot = p.sum(0)
        prob[ids] = np.where(mask, p/np.maximum(tot, 1e-9), 0)
    label = np.where(rcord, prob.argmax(0), -1)

    gm_before = dice(pg > 0.5, rg)
    gm_after = dice(wg > 0.5, rg)
    tracts = []
    for i, name, _ in names:
        side, short = name.split()[1], " ".join(name.split()[2:])
        region = label == i
        tracts.append({"id": i, "name": name, "short": short, "side": side,
                       "tissue": "WM" if i < 30 else "GM", "colour": COLOURS[short],
                       "area_argmax_mm2": round(float(region.sum()*H*H), 4),
                       "area_prob_mm2": round(float(prob[i].sum()*H*H), 4),
                       "polygons_uv_mm": polygons(region)})
    rep = {"pam50_level": args.level, "pam50_slice": k, "rado_section_s_mm": args.s,
           "rado_section_point_mm": c0.round(4).tolist(),
           "frame": "u = +x (patient left), v = posterior, perpendicular to FIBRE_CENTRELINE; "
                    "(u, v) are offsets from the centreline point at any arc length",
           "centreline": os.path.relpath(C.FIBRE_CENTRELINE, ROOT),
           "grid_mm": H, "sct_param": os.environ.get("SCS_SCT_PARAM", PARAM), "blur_mm": float(os.environ.get("SCS_BLUR_MM", BLUR_MM)),
           "grey_dice_before": round(gm_before, 4), "grey_dice_after": round(gm_after, 4),
           "cord_area_mm2": round(float(rcord.sum()*H*H), 3),
           "tracts": tracts, "method": __doc__}
    with open(os.path.join(out, "tracts.json"), "w") as f:
        json.dump(rep, f, indent=1)
    np.savez_compressed(os.path.join(out, "tract_map.npz"), u=g, v=g, prob=prob.astype(np.float32),
                        label=label.astype(np.int8), white=rw, grey=rg)

    # figure: PAM50 source labels (unwarped), RADO target, warped labels
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap([COLOURS[t["short"]] for t in tracts])
    src_lab = np.where((pc > 0.5), maps.argmax(0), -1)
    fig, axs = plt.subplots(1, 3, figsize=(15, 5.4))
    for ax, lab, title in ((axs[0], src_lab, "PAM50 %s (slice %d), centred" % (args.level, k)),
                           (axs[2], label, "atlas warped onto RADO (grey Dice %.2f -> %.2f)" % (gm_before, gm_after))):
        ax.imshow(np.ma.masked_less(lab, 0).T, cmap=cmap, vmin=-0.5, vmax=35.5, extent=[g[0], g[-1], g[-1], g[0]],
                  interpolation="nearest")
        ax.contour(g, g, rg.T.astype(float), [0.5], colors="k", linewidths=0.8)
        ax.contour(g, g, rcord.T.astype(float), [0.5], colors="k", linewidths=0.8)
        ax.set_title(title, fontsize=10)
    axs[1].imshow(dst.T, cmap="gray_r", extent=[g[0], g[-1], g[-1], g[0]])
    axs[1].set_title("RADO section (s = %.0f mm on the extended centreline)" % args.s, fontsize=10)
    for ax in axs:
        ax.set_xlabel("u, mm (+ = patient left)"); ax.set_ylabel("v, mm (+ = posterior, down)")
        ax.set_aspect("equal")
    seen = {}
    for t in tracts:
        seen.setdefault(t["short"], t["colour"])
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in seen.values()]
    fig.legend(handles, list(seen), loc="lower center", ncol=6, fontsize=7)
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.savefig(os.path.join(out, "tract_map.png"), dpi=140)

    print("grey Dice %.3f -> %.3f; cord %.1f mm2; written to %s" % (gm_before, gm_after, rep["cord_area_mm2"], out))
    for t in tracts:
        if t["side"] == "left":
            print("  %-60s %6.2f mm2" % (t["short"], t["area_argmax_mm2"]))


if __name__ == "__main__":
    main()
