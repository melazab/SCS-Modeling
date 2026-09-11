"""Sagittal and axial slices of the bipolar 1 A voltage field through the lead."""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
from field import TetField

# Outline colours for the tissue boundaries drawn over the field. Chosen for
# contrast against a red/blue diverging map, not to match the FreeCAD document.
TISSUE_COLOR = {"grey": "#00a0a0", "white": "#00d2d2", "csf": "#7b3fbf",
                "dura": "#111111", "epidural": "#9a7b00", "insulator": "#00ff90"}
TISSUE_LABEL = {"grey": "grey matter", "white": "white matter", "csf": "CSF",
                "dura": "dura (0.5 mm)", "epidural": "epidural fat",
                "insulator": "lead body"}


def slab(f, lab_vals, fixed_axis, fixed, u_rng, v_rng, n=420):
    ua = np.linspace(*u_rng, n)
    va = np.linspace(*v_rng, int(n * (v_rng[1] - v_rng[0]) / (u_rng[1] - u_rng[0])))
    U, Vv = np.meshgrid(ua, va, indexing="ij")
    P = np.empty((U.size, 3))
    ax = [a for a in range(3) if a != fixed_axis]
    P[:, ax[0]] = U.ravel()
    P[:, ax[1]] = Vv.ravel()
    P[:, fixed_axis] = fixed
    ti, w = f.locate(P)
    val = np.full(U.size, np.nan)
    ok = ti >= 0
    vn = f.values[0][f.tets[ti[ok]]]
    val[ok] = np.einsum("nk,nk->n", vn, w[ok])
    tis = np.where(ok, lab_vals[np.clip(ti, 0, None)], -1)
    return ua, va, val.reshape(U.shape), tis.reshape(U.shape)


def main():
    d = np.load(os.path.join(C.OUT, "solution.npz"), allow_pickle=True)
    order = list(d["order"])
    f = TetField(d["nodes"], d["tets"], d["V"].astype(np.float64)[None, :])
    lab = d["label"]

    c3 = np.array([55.60, 87.49, 132.83])

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 8.4),
                             gridspec_kw=dict(width_ratios=[1.55, 1]))

    # --- sagittal, x = lead axis ---
    ua, va, val, tis = slab(f, lab, 0, 55.60, (70.0, 92.0), (112.0, 166.0), 300)
    amax = np.nanpercentile(np.abs(val), 99.8)
    norm = SymLogNorm(linthresh=max(amax / 300, 1e-4), vmin=-amax, vmax=amax, base=10)
    # log-spaced iso-potential lines over the range actually present
    dec = amax * 10.0 ** -np.arange(0, 3.5, 0.5)
    lvl = np.sort(np.concatenate([-dec, dec]))
    ax = axes[0]
    im = ax.pcolormesh(ua, va, val.T, cmap="RdBu_r", norm=norm, shading="auto")
    ax.contour(ua, va, val.T, levels=np.sort(lvl), colors="k", linewidths=0.35, alpha=.5)
    for name, col in TISSUE_COLOR.items():
        if name not in order:
            continue
        ax.contour(ua, va, (tis == order.index(name)).T.astype(float),
                   levels=[0.5], colors=[col], linewidths=1.1)
    ax.set_xlabel("y  (mm, ventral → dorsal)")
    ax.set_ylabel("z  (mm, rostral → up)")
    ax.set_title("Sagittal slice, x = 55.6 mm (lead axis)")
    ax.set_aspect("equal")
    ax.set_xlim(70.0, 98.0)   # headroom on the right for the contact labels
    for i, (z0, z1) in enumerate([(123.25, 126.41), (127.26, 130.40), (131.27, 134.39),
                                  (135.28, 138.38), (139.29, 142.37), (143.30, 146.36),
                                  (147.31, 150.35), (151.32, 154.34)], start=1):
        tag = {3: "+1 A", 5: "−1 A"}.get(i, "float")
        ax.text(92.6, (z0 + z1) / 2, "c%d  %s" % (i, tag), fontsize=8,
                va="center", ha="left", clip_on=False,
                color="firebrick" if i == 3 else ("navy" if i == 5 else "0.35"))

    # --- axial, z = centre of contact 3 ---
    ua2, va2, val2, tis2 = slab(f, lab, 2, c3[2], (44.0, 69.0), (71.0, 91.0), 460)
    ax = axes[1]
    ax.pcolormesh(ua2, va2, val2.T, cmap="RdBu_r", norm=norm, shading="auto")
    ax.contour(ua2, va2, val2.T, levels=np.sort(lvl), colors="k", linewidths=0.35, alpha=.5)
    for name, col in TISSUE_COLOR.items():
        if name not in order:
            continue
        ax.contour(ua2, va2, (tis2 == order.index(name)).T.astype(float),
                   levels=[0.5], colors=[col], linewidths=1.1)
    ax.set_xlabel("x  (mm, mediolateral)")
    ax.set_ylabel("y  (mm, dorsal → up)")
    ax.set_title("Axial slice, z = %.1f mm (centre of contact 3)" % c3[2])
    ax.set_aspect("equal")

    from matplotlib.lines import Line2D
    axes[1].legend(handles=[Line2D([], [], color=TISSUE_COLOR[k], lw=1.6,
                                   label=TISSUE_LABEL[k])
                            for k in ("epidural", "dura", "csf", "white", "grey",
                                      "insulator")],
                   loc="upper center", bbox_to_anchor=(0.5, -0.16), fontsize=8,
                   ncol=3, framealpha=0.95, title="tissue boundaries",
                   title_fontsize=8)
    cb = fig.colorbar(im, ax=axes, fraction=0.03, pad=0.02)
    cb.set_label("potential  V  (volts, bipolar 1 A: +1 A c3, −1 A c5)")
    fig.suptitle("RADO-SCS stripped model — gmsh + P1 FEM, div(σ∇V)=0, "
                 "bipolar 1 A dorsal epidural", fontsize=11)
    out = os.path.join(C.OUT, "voltage_slices.png")
    fig.savefig(out, dpi=155, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
