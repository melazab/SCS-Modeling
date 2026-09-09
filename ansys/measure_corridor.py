#!/usr/bin/env python3
"""Measure the dorsal and ventral epidural corridors as a function of height.

The T8-T10 spine in this model is kyphotic, so the epidural space does not sit
at a constant anterior-posterior position: near midline the dorsal channel
centre moves from y ~= 81.5 mm at z = 95 to y ~= 87.5 mm at z = 135, about 6 mm
over the length of a clinical 8-contact lead. A straight lead at fixed y
therefore walks out of the epidural space -- which is presumably why RADO's own
lead is curved.

This measures the corridor from the STL surfaces and fits a smooth centreline
that ansys/make_scs_lead.py can sweep a lead along, so the generated lead
follows the canal instead of cutting through the dura.

Frame (see ansys/check_laterality.py): +X anatomical left, midline x = 56.60;
+Y posterior/dorsal; +Z rostral. Distances in mm.

    python3 ansys/measure_corridor.py [--half-width 2.0] [--degree 3]

Writes ansys/epidural_corridor.json.
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from find_floating_bodies import read_stl_vertices  # noqa: E402

MIDLINE_X = 56.60


def band_profile(epi, dura, x_half, z_lo, z_hi, step):
    """Per z-slice: (z, ventral centre, ventral width, dorsal centre, dorsal width)."""
    rows = []
    z = z_lo
    while z < z_hi:
        e = epi[(np.abs(epi[:, 0] - MIDLINE_X) < x_half) &
                (epi[:, 2] >= z) & (epi[:, 2] < z + step)]
        d = dura[(np.abs(dura[:, 0] - MIDLINE_X) < x_half) &
                 (dura[:, 2] >= z) & (dura[:, 2] < z + step)]
        if len(e) >= 20 and len(d) >= 20:
            e_lo, e_hi = float(e[:, 1].min()), float(e[:, 1].max())
            d_lo, d_hi = float(d[:, 1].min()), float(d[:, 1].max())
            rows.append(dict(z=z + 0.5 * step,
                             ventral_centre=0.5 * (e_lo + d_lo),
                             ventral_width=d_lo - e_lo,
                             dorsal_centre=0.5 * (d_hi + e_hi),
                             dorsal_width=e_hi - d_hi))
        z += step
    return rows


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument("--stl-dir", default=os.path.join(here, "..", "STL_files"))
    ap.add_argument("--half-width", type=float, default=2.0,
                    help="how far either side of midline to sample, mm")
    ap.add_argument("--step", type=float, default=2.0, help="z slice thickness, mm")
    ap.add_argument("--degree", type=int, default=3, help="polynomial degree for the fit")
    ap.add_argument("--out", default=os.path.join(here, "epidural_corridor.json"))
    args = ap.parse_args()

    epi = read_stl_vertices(os.path.join(args.stl_dir, "T8-10 - neuro_EpiduralSpace-1.STL"))
    dura = read_stl_vertices(os.path.join(args.stl_dir, "T8-10 - neuro_Meninges-1.STL"))
    z_lo = max(epi[:, 2].min(), dura[:, 2].min())
    z_hi = min(epi[:, 2].max(), dura[:, 2].max())
    rows = band_profile(epi, dura, args.half_width, z_lo, z_hi, args.step)
    if len(rows) < args.degree + 2:
        raise SystemExit("too few usable z slices (%d)" % len(rows))

    z = np.array([r["z"] for r in rows])
    out = {"frame": {"midline_x": MIDLINE_X, "x_is": "anatomical left (+)",
                     "y_is": "posterior/dorsal (+)", "z_is": "rostral (+)"},
           "z_range": [float(z.min()), float(z.max())],
           "sample_step_mm": args.step, "fit_degree": args.degree, "sides": {}}

    print("%-8s %-10s %-9s %-10s %-9s" %
          ("z", "vent ctr", "vent w", "dors ctr", "dors w"))
    for r in rows:
        print("%-8.1f %-10.2f %-9.2f %-10.2f %-9.2f"
              % (r["z"], r["ventral_centre"], r["ventral_width"],
                 r["dorsal_centre"], r["dorsal_width"]))

    for side in ("dorsal", "ventral"):
        c = np.array([r["%s_centre" % side] for r in rows])
        w = np.array([r["%s_width" % side] for r in rows])
        keep = w > 0.2                     # slices where a channel actually exists
        coef = np.polyfit(z[keep], c[keep], args.degree)
        resid = c[keep] - np.polyval(coef, z[keep])
        out["sides"][side] = {
            "centre_poly_coef_high_to_low": [float(v) for v in coef],
            "fit_rms_mm": float(np.sqrt((resid ** 2).mean())),
            "fit_max_abs_mm": float(np.abs(resid).max()),
            "width_min_mm": float(w[keep].min()),
            "width_median_mm": float(np.median(w[keep])),
            "usable_z": [float(z[keep].min()), float(z[keep].max())],
            "n_slices": int(keep.sum()),
        }
        s = out["sides"][side]
        print("\n%s: centreline fit deg %d, rms %.3f mm, max %.3f mm over z %.1f..%.1f"
              % (side, args.degree, s["fit_rms_mm"], s["fit_max_abs_mm"],
                 s["usable_z"][0], s["usable_z"][1]))
        print("   channel width: min %.2f mm, median %.2f mm  -> widest lead that fits: %.2f mm"
              % (s["width_min_mm"], s["width_median_mm"], s["width_min_mm"]))

    json.dump(out, open(args.out, "w"), indent=2)
    print("\nwrote %s" % args.out)


if __name__ == "__main__":
    main()
