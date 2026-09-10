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
curved. The fitted polynomials here are what src/freecad/make_scs_lead.py sweeps along.

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
import os
import shlex
import sys

import FreeCAD as App
import Mesh

MIDLINE_X = 56.60
HERE = os.path.dirname(os.path.abspath(__file__))
STL_DIR = os.path.join(HERE, "..", "..", "STL_files")


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

    epi = Mesh.Mesh(os.path.join(STL_DIR, "T8-10 - neuro_EpiduralSpace-1.STL"))
    bb = epi.BoundBox
    lines = ["epidural corridor at x = %.2f (midline), ray cast along +Y" % MIDLINE_X,
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

    zs = [r["z"] for r in rows]
    out = {"frame": {"midline_x": MIDLINE_X, "x_is": "anatomical left (+)",
                     "y_is": "posterior/dorsal (+)", "z_is": "rostral (+)"},
           "method": "ray cast along +Y through the epidural mesh at x=midline",
           "sample_step_mm": step, "fit_degree": degree,
           "z_range": [min(zs), max(zs)], "n_slices": len(rows), "sides": {}}

    for side in ("dorsal", "ventral"):
        cs = [r["%s_centre" % side] for r in rows]
        ts = [r["%s_thickness" % side] for r in rows]
        coef = polyfit(zs, cs, degree)
        resid = [c - polyval(coef, z) for z, c in zip(zs, cs)]
        rms = (sum(r * r for r in resid) / len(resid)) ** 0.5
        out["sides"][side] = {
            "centre_poly_coef_high_to_low": coef,
            "fit_rms_mm": rms,
            "fit_max_abs_mm": max(abs(r) for r in resid),
            "thickness_min_mm": min(ts),
            "thickness_median_mm": sorted(ts)[len(ts) // 2],
            "usable_z": [min(zs), max(zs)],
        }
        s = out["sides"][side]
        lines += ["", "%s: centreline y(z) deg %d, rms %.3f mm, max %.3f mm"
                  % (side, degree, s["fit_rms_mm"], s["fit_max_abs_mm"]),
                  "   fat thickness min %.2f mm, median %.2f mm"
                  % (s["thickness_min_mm"], s["thickness_median_mm"]),
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
