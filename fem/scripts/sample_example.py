"""Worked example: sample the FEM potential along a dorsal-column fibre.

This is the handover point to the NEURON stage. It shows both routes:

  A. trilinear interpolation of fem/out/field_grid.npz  -- fast, what
     docs/neuron_plan.md chunk 3 asks for
  B. direct barycentric sampling of the tetrahedral mesh via field.TetField --
     exact at arbitrary points, no grid resampling error

and it prints the activating function, d2V/ds2 along the fibre, because that --
not V itself -- is what drives an axon and what will show up any interpolation
noise first (neuron_plan's "if it is noisy, the FEM grid is too coarse").

Units here stay mm and volts. The NEURON side converts once, to um and mV.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
from field import TetField


def trilinear(g, pts):
    origin, sp, shape = g["origin"], float(g["spacing"]), g["shape"]
    phi = g["phi"]
    t = (pts - origin) / sp
    i0 = np.floor(t).astype(int)
    f = t - i0
    out = np.zeros((phi.shape[0], len(pts)))
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                idx = i0 + [dx, dy, dz]
                ok = np.all((idx >= 0) & (idx < shape), axis=1)
                w = (np.where(dx, f[:, 0], 1 - f[:, 0])
                     * np.where(dy, f[:, 1], 1 - f[:, 1])
                     * np.where(dz, f[:, 2], 1 - f[:, 2]))
                ii = np.clip(idx, 0, np.asarray(shape) - 1)
                out[:, ok] += w[ok] * phi[:, ii[ok, 0], ii[ok, 1], ii[ok, 2]]
    return out


def main():
    # A straight rostrocaudal fibre in dorsal white matter, 1 mm deep to the
    # dorsal cord surface and on the lead's side of the midline. Chunk 4 will
    # replace this with real trajectories from the anatomy.
    z = np.arange(112.0, 156.0, 0.5)
    pts = np.stack([np.full_like(z, 55.6), np.full_like(z, 82.3), z], axis=1)

    g = np.load(os.path.join(C.OUT, "field_grid.npz"), allow_pickle=True)
    phi_grid = trilinear(g, pts)
    Vgrid = phi_grid[C.SOURCE_CONTACT - 1] - phi_grid[C.SINK_CONTACT - 1]

    d = np.load(os.path.join(C.OUT, "solution.npz"), allow_pickle=True)
    f = TetField(d["nodes"], d["tets"], d["phi"].astype(np.float64))
    phi_tet = f(pts)
    Vtet = phi_tet[C.SOURCE_CONTACT - 1] - phi_tet[C.SINK_CONTACT - 1]

    ds = z[1] - z[0]
    act = np.full_like(Vtet, np.nan)
    act[1:-1] = (Vtet[2:] - 2 * Vtet[1:-1] + Vtet[:-2]) / (ds * 1e-3) ** 2

    print("dorsal-column fibre at x=55.6, y=82.3 mm, bipolar 1 A (+c3 / -c5)")
    print("%8s %12s %12s %12s %14s" % ("z mm", "V grid (V)", "V tets (V)",
                                       "diff (V)", "d2V/ds2 V/m2"))
    for i in range(0, len(z), 4):
        print("%8.1f %12.4f %12.4f %12.2e %14.4g"
              % (z[i], Vgrid[i], Vtet[i], Vgrid[i] - Vtet[i], act[i]))
    dif = np.abs(Vgrid - Vtet)
    print("\ngrid vs tet sampling: max %.3e V, rms %.3e V, over a %.3f V swing"
          % (np.nanmax(dif), np.sqrt(np.nanmean(dif ** 2)), np.ptp(Vtet)))
    k = np.nanargmax(np.abs(act))
    print("activating function peaks at z = %.1f mm, %.4g V/m^2" % (z[k], act[k]))


if __name__ == "__main__":
    main()
