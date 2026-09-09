#!/usr/bin/env python3
"""Plot a voltage-spread figure from a MAPDL slice export.

Input is the CSV written by the `*VWRITE` block in the decks
(`ansys/testA_mapdl_cylinder/testA.dat`, `ansys/bigdeck_recondition/`):
three columns, no header -- in-plane coordinate 1, in-plane coordinate 2, VOLT.

Uses matplotlib's own Delaunay triangulation (matplotlib.tri), so it needs only
numpy + matplotlib -- scipy is NOT installed on this machine.

    python3 plot_voltage_slice.py slice.csv -o fig.png --units mm --title "..."
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import TwoSlopeNorm, SymLogNorm

UNIT_SCALE = {"m": 1.0, "mm": 1e-3, "um": 1e-6}  # CSV unit -> metres


def load(path, unit):
    d = np.loadtxt(path, delimiter=",")
    if d.ndim != 2 or d.shape[1] < 3:
        raise SystemExit("expected 3 columns (c1, c2, volt), got shape %s" % (d.shape,))
    a, b, v = d[:, 0], d[:, 1], d[:, 2]
    # collapse duplicate coordinates (mid-side nodes can coincide in projection)
    scale = UNIT_SCALE[unit] / UNIT_SCALE["mm"]  # render in mm
    return a * scale, b * scale, v


def main():
    p = argparse.ArgumentParser()
    p.add_argument("csv")
    p.add_argument("-o", "--out", default="voltage_slice.png")
    p.add_argument("--units", default="m", choices=sorted(UNIT_SCALE),
                   help="units of the coordinates in the CSV (default m)")
    p.add_argument("--title", default="Bipolar SCS: voltage spread")
    p.add_argument("--xlabel", default="x [mm]")
    p.add_argument("--ylabel", default="z (rostro-caudal) [mm]")
    p.add_argument("--clip", type=float, default=99.5,
                   help="percentile of |V| used for the colour limits (default 99.5)")
    args = p.parse_args()

    x, y, v = load(args.csv, args.units)
    print("%d points | x %.2f..%.2f mm | y %.2f..%.2f mm | V %.4g..%.4g"
          % (len(v), x.min(), x.max(), y.min(), y.max(), v.min(), v.max()))

    tri = mtri.Triangulation(x, y)
    # drop slivers produced by triangulating a projected point cloud
    xt, yt = x[tri.triangles], y[tri.triangles]
    area = 0.5 * np.abs((xt[:, 1] - xt[:, 0]) * (yt[:, 2] - yt[:, 0])
                        - (xt[:, 2] - xt[:, 0]) * (yt[:, 1] - yt[:, 0]))
    edge = np.max(np.hypot(np.diff(np.c_[xt, xt[:, :1]], axis=1),
                           np.diff(np.c_[yt, yt[:, :1]], axis=1)), axis=1)
    tri.set_mask((area < 1e-9) | (edge > np.percentile(edge, 99.0) * 3))

    lim = np.percentile(np.abs(v), args.clip)
    if lim <= 0:
        lim = max(abs(v).max(), 1e-12)

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.6), constrained_layout=True)

    ax = axes[0]
    norm = TwoSlopeNorm(vmin=-lim, vcenter=0.0, vmax=lim)
    cf = ax.tricontourf(tri, v, levels=np.linspace(-lim, lim, 41),
                        cmap="RdBu_r", norm=norm, extend="both")
    ax.tricontour(tri, v, levels=np.linspace(-lim, lim, 21),
                  colors="k", linewidths=0.25, alpha=0.35)
    fig.colorbar(cf, ax=ax, label="V [volts]")
    ax.set_title("linear scale")

    ax = axes[1]
    small = lim / 1e3
    cf2 = ax.tricontourf(tri, v, levels=41, cmap="RdBu_r",
                         norm=SymLogNorm(linthresh=small, vmin=-lim, vmax=lim, base=10),
                         extend="both")
    fig.colorbar(cf2, ax=ax, label="V [volts], symlog")
    ax.set_title("symmetric-log scale (shows far-field spread)")

    for ax in axes:
        ax.set_aspect("equal")
        ax.set_xlabel(args.xlabel)
        ax.set_ylabel(args.ylabel)

    fig.suptitle(args.title)
    fig.savefig(args.out, dpi=160)
    print("wrote %s" % args.out)


if __name__ == "__main__":
    main()
