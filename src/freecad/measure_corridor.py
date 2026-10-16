"""Measure the dorsal and ventral epidural corridors by ray casting.

Run with FreeCAD (stdout is swallowed by freecadcmd, so results also go to a
report file and to src/freecad/epidural_corridor.json):

    freecadcmd src/freecad/measure_corridor.py
    MEASURE_CORRIDOR_ARGS="--step 2" freecadcmd src/freecad/measure_corridor.py

Why ray casting and not surface arithmetic
------------------------------------------
RADO's compartments are HOLLOW SHELLS that tile space -- each body has the
inner ones carved out of it. Verified by ray-casting a point at the centre of a
cord cross-section: it registers 3 forward crossings of the white-matter surface
(odd => inside the material) but exactly 2 for CSF, dura and the epidural body
(even => in their central cavity, i.e. outside their material).

An earlier version of this script measured the corridor as
`epidural.ymax - dura.ymax`, i.e. the distance between two OUTER surfaces, and
concluded there was only ~0.03 mm of epidural fat. That was wrong -- it is not
the thickness of anything. Casting a ray along +Y and reading off consecutive
entry/exit pairs gives the true material intervals:

    z = 95 mm, x = midline:  ventral fat y 67.60..69.34 (1.74 mm)
                             dorsal  fat y 80.12..82.47 (2.35 mm)

which is anatomically sensible -- the dorsal epidural space is the thicker of
the two -- and means a 1.3 mm clinical lead fits on either side.

The canal also migrates posteriorly with height (kyphosis): the dorsal centre
runs y ~= 81.3 at z = 95 to y ~= 88.4 at z = 145, about 7 mm over the length of
an 8-contact lead. A lead has to follow that, which is why RADO's own lead is
curved. The stations measured here (every 2 mm, written to epidural_corridor.json)
are what a lead's local cubic is fitted to; see main() for why the extended
T6-T12 canal is stored as stations rather than one polynomial, and why the
thickness is corrected to the canal normal. The anatomy set is anatomy.py's.

Frame (see src/freecad/check_laterality.py): +X anatomical left, midline x = 56.60;
+Y posterior/dorsal; +Z rostral. Millimetres throughout.

IMPORTING THIS MODULE IS SAFE
-----------------------------
src/freecad/build_lead_config.py imports intervals(), polyfit() and polyval()
from here, so that ray casting the epidural mesh and fitting a centreline have
one implementation each. Until the guard at the bottom existed, importing this
file re-ran the whole measurement and REWROTE epidural_corridor.json as a side
effect -- freecadcmd imports a script rather than exec'ing it, so the usual
__main__ guard never fires and main() had to be called at module level.
"""
import json
import math
import os
import shlex
import sys

import FreeCAD as App
import Mesh

MIDLINE_X = 56.60
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import anatomy  # noqa: E402

STL_DIR = anatomy.STL_DIR


def polyfit(xs, ys, deg):
    """Least-squares polynomial fit via normal equations (no numpy dependency
    inside freecadcmd). Returns coefficients highest power first."""
    n = deg + 1
    A = [[sum(x ** (i + j) for x in xs) for j in range(n)] for i in range(n)]
    b = [sum(y * x ** i for x, y in zip(xs, ys)) for i in range(n)]
    for c in range(n):                       # Gaussian elimination w/ partial pivot
        p = max(range(c, n), key=lambda r: abs(A[r][c]))
        A[c], A[p] = A[p], A[c]
        b[c], b[p] = b[p], b[c]
        for r in range(c + 1, n):
            f = A[r][c] / A[c][c]
            for k in range(c, n):
                A[r][k] -= f * A[c][k]
            b[r] -= f * b[c]
    coef = [0.0] * n
    for r in range(n - 1, -1, -1):
        s = b[r] - sum(A[r][k] * coef[k] for k in range(r + 1, n))
        coef[r] = s / A[r][r]
    return list(reversed(coef))              # highest power first


def polyval(coef, x):
    v = 0.0
    for c in coef:
        v = v * x + c
    return v


def local_fit(stations, z_lo, z_hi, degree=3):
    """y(z) through the stations with z_lo <= z <= z_hi, as a callable.

    Fitted in z - z_ref, z_ref the span's middle: at z ~ 280 the raw powers
    1, z, z^2, z^3 are nearly collinear and the normal equations lose every
    significant digit. Raises if the span holds too few stations."""
    pts = [(st[0], st[1]) for st in stations if z_lo <= st[0] <= z_hi]
    if len(pts) < degree + 2:
        raise ValueError("only %d measured stations in z %.1f..%.1f" % (len(pts), z_lo, z_hi))
    z_ref = 0.5 * (z_lo + z_hi)
    coef = polyfit([z - z_ref for z, _ in pts], [y for _, y in pts], degree)
    return lambda z: polyval(coef, z - z_ref)


def slope_at(stations, z, half_window=6.0):
    """dy/dz of the channel centre at z, from a local quadratic."""
    pts = [(st[0], st[1]) for st in stations if abs(st[0] - z) <= half_window]
    if len(pts) < 4:
        return 0.0
    coef = polyfit([q - z for q, _ in pts], [y for _, y in pts], 2)
    return coef[1]


def accept_single(segs, y_pred, t_ref, tol_mm=1.0, tol_frac=0.15):
    """A slice that crosses only ONE channel, accepted as the channel whose
    curve it continues, else None.

    RADO's canal ends in a plane tilted ~42 deg (extended model) or ~18 deg
    (RADO's own top), so near the top the ventral fat has already ended while
    the dorsal fat runs on for another ~10 mm. A +Y ray there gives one
    interval, not two. It is that channel if its centre is within tol_mm of
    where the channel's own fitted curve says it should be AND its thickness
    is within tol_frac of the channel's; anything else (the cap cutting the
    fat obliquely, a root sleeve) is refused."""
    if len(segs) != 1:
        return None
    lo, hi, t = segs[0]
    if abs(0.5 * (lo + hi) - y_pred) <= tol_mm and abs(t - t_ref) <= tol_frac * t_ref:
        return segs[0]
    return None


def extend_single(mesh, x, stations, step, z_stop, direction):
    """Continue one channel's stations [(z, y, t)] past its last two-channel
    slice, one step at a time, while accept_single() agrees. Returns the
    added stations, nearest first."""
    added = []
    ordered = sorted(stations)
    edge = ordered[-1] if direction > 0 else ordered[0]
    z = edge[0] + direction * step
    while (z - z_stop) * direction <= 0:
        near = sorted(stations + added, key=lambda st: abs(st[0] - z))[:15]
        f = local_fit(near, min(st[0] for st in near), max(st[0] for st in near) + 1e-9, 2)
        t_ref = sorted(st[2] for st in near)[len(near) // 2]
        hit = accept_single(intervals(mesh, x, z), f(z), t_ref)
        if hit is None:
            break
        added.append((z, 0.5 * (hit[0] + hit[1]), hit[2]))
        z += direction * step
    return added


def intervals(mesh, x, z):
    """Material intervals along +Y through (x, ., z), as (y_in, y_out, thickness)."""
    hits = mesh.foraminate((x, 0.0, z), (0.0, 1.0, 0.0))
    ys = sorted(v[1] for v in hits.values())
    return [(ys[i], ys[i + 1], ys[i + 1] - ys[i]) for i in range(0, len(ys) - 1, 2)]


def main():
    argv = shlex.split(os.environ.get("MEASURE_CORRIDOR_ARGS", ""))
    step = 2.0
    degree = 3
    if "--step" in argv:
        step = float(argv[argv.index("--step") + 1])
    if "--degree" in argv:
        degree = int(argv[argv.index("--degree") + 1])

    epi_name = anatomy.stl_name("neuro_EpiduralSpace-1")
    epi = Mesh.Mesh(os.path.join(STL_DIR, epi_name))
    bb = epi.BoundBox
    lines = ["epidural corridor at x = %.2f (midline), ray cast along +Y" % MIDLINE_X,
             "anatomy: %s (%s)" % (anatomy.NAME, epi_name),
             "", "%-8s %-24s %-24s" % ("z", "ventral fat", "dorsal fat")]

    rows = []
    z = bb.ZMin + step
    while z < bb.ZMax - step:
        segs = intervals(epi, MIDLINE_X, z)
        # a clean slice gives exactly two intervals: ventral then dorsal
        if len(segs) == 2:
            (va, vb, vt), (da, db, dt) = segs
            rows.append(dict(z=z,
                             ventral_centre=0.5 * (va + vb), ventral_thickness=vt,
                             dorsal_centre=0.5 * (da + db), dorsal_thickness=dt))
            lines.append("%-8.1f %6.2f..%6.2f (%4.2f)   %6.2f..%6.2f (%4.2f)"
                         % (z, va, vb, vt, da, db, dt))
        else:
            lines.append("%-8.1f  %d intervals - skipped (branching or outside the sac)"
                         % (z, len(segs)))
        z += step

    if len(rows) < degree + 2:
        lines.append("\nToo few clean slices (%d) to fit." % len(rows))
        open(os.path.join(HERE, "corridor_report.txt"), "w").write("\n".join(lines))
        return

    # THE STATIONS ARE THE RESULT. RADO's own canal (z 61-159) was one cubic
    # to 0.02 mm, and this file used to hold only that cubic. The extended
    # canal turns ~43 deg between T12 and T6, which no single low-order
    # polynomial follows, so the stations themselves are stored and every
    # lead fits its own cubic to the ones along it (local_fit).
    #
    # THICKNESS PERPENDICULAR TO THE CANAL. A +Y ray crosses a channel tilted
    # by theta from z obliquely, reading t / cos(theta); the 20-23 deg tilts at
    # the two ends of the extended canal would overstate it by 6-9 %. Each
    # station's thickness is corrected with the local slope of its own centre
    # line. (The ray's MIDPOINT is unaffected: it lies on the channel's
    # mid-surface whatever the angle.)
    zs = [r["z"] for r in rows]
    out = {"frame": {"midline_x": MIDLINE_X, "x_is": "anatomical left (+)",
                     "y_is": "posterior/dorsal (+)", "z_is": "rostral (+)"},
           "anatomy": anatomy.NAME, "epidural_stl": epi_name,
           "method": "ray cast along +Y through the epidural mesh at x=midline; "
                     "thickness corrected to the canal normal with the local slope",
           "station_columns": ["z", "y_centre", "thickness_along_y", "thickness_normal"],
           "sample_step_mm": step, "local_fit_degree": degree,
           "z_range": [min(zs), max(zs)], "n_slices": len(rows), "sides": {}}

    for side in ("dorsal", "ventral"):
        raw = [(r["z"], r["%s_centre" % side], r["%s_thickness" % side]) for r in rows]
        # past the ends of the two-channel run, where only this channel remains
        up = extend_single(epi, MIDLINE_X, raw, step, bb.ZMax, +1)
        down = extend_single(epi, MIDLINE_X, raw, step, bb.ZMin, -1)
        raw = sorted(raw + up + down)
        if up or down:
            lines.append("%s: %d single-channel slices accepted above z %.1f, %d below z %.1f"
                         % (side, len(up), max(zs), len(down), min(zs)))
        stations = []
        worst_tilt = 0.0
        for z, c, t in raw:
            k = slope_at(raw, z)
            worst_tilt = max(worst_tilt, abs(math.degrees(math.atan(k))))
            stations.append([z, c, t, t * math.cos(math.atan(k))])
        # quality of the local fits a lead will use: a cubic over each 60 mm window
        worst = 0.0
        for z, c, _t, _n in stations:
            try:
                f = local_fit(stations, z - 30.0, z + 30.0, degree)
            except ValueError:
                continue
            worst = max(worst, abs(f(z) - c))
        tn = [st[3] for st in stations]
        out["sides"][side] = {
            "stations": stations,
            "local_fit_max_abs_mm": worst,
            "thickness_min_mm": min(tn),
            "thickness_median_mm": sorted(tn)[len(tn) // 2],
            "thickness_min_at_z": stations[tn.index(min(tn))][0],
            "max_tilt_deg": worst_tilt,
            "usable_z": [stations[0][0], stations[-1][0]],
        }
        s = out["sides"][side]
        lines += ["", "%s: %d stations, cubic over any 60 mm window fits to %.3f mm; "
                  "canal tilts up to %.1f deg from z"
                  % (side, len(stations), s["local_fit_max_abs_mm"], s["max_tilt_deg"]),
                  "   fat thickness normal to the canal: min %.2f mm (z %.0f), median %.2f mm"
                  % (s["thickness_min_mm"], s["thickness_min_at_z"], s["thickness_median_mm"]),
                  "   widest lead that fits: %.2f mm" % s["thickness_min_mm"]]

    json.dump(out, open(os.path.join(HERE, "epidural_corridor.json"), "w"), indent=2)
    lines += ["", "wrote src/freecad/epidural_corridor.json"]
    open(os.path.join(HERE, "corridor_report.txt"), "w").write("\n".join(lines) + "\n")
    sys.stdout.write("\n".join(lines) + "\n")


def run_as_freecadcmd_script():
    """True when a FreeCAD interpreter was handed THIS file to run.

    Same guard as apply_colors.py. Without it, `from measure_corridor import
    intervals` would re-measure the corridor and overwrite
    epidural_corridor.json -- see the module docstring.
    """
    if len(sys.argv) < 2:
        return False
    if not os.path.basename(sys.argv[0]).lower().startswith("freecad"):
        return False
    return os.path.abspath(sys.argv[1]) == os.path.abspath(__file__)


if __name__ == "__main__":
    main()
elif run_as_freecadcmd_script():
    main()
