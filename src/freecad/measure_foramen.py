"""Measure the eight foraminal corridors -- the paths a DRG lead can take.

Run with FreeCAD (stdout is swallowed by freecadcmd, so results also go to a
report file and to src/freecad/foraminal_corridors.json):

    freecadcmd src/freecad/measure_foramen.py
    MEASURE_FORAMEN_ARGS="--step 0.5 --degree 4" freecadcmd src/freecad/measure_foramen.py

WHY A DRG LEAD NEEDS ITS OWN MEASUREMENT
----------------------------------------
A dorsal or ventral epidural lead runs rostro-caudally along the canal, and
src/freecad/measure_corridor.py measures that canal as a channel y(z) at fixed
x. A percutaneous DRG lead does not go there at all. It leaves the canal at one
level and runs LATERALLY, out through the neuroforamen, so that its contact
array ends up beside a specific dorsal root ganglion. Its long axis is x, not z,
and the two coordinates that have to follow a curve are y and z. Sweeping such a
lead along the canal centreline would put it in the wrong plane entirely.

So this is the DRG analogue of measure_corridor.py: same ray casting, same
polynomial fit, different axis -- a centreline (y(u), z(u)) parametrised by u,
the LATERAL DISTANCE FROM THE MIDLINE, for each of the eight foramina.

WHAT IS BEING RAY CAST
----------------------
RADO models the tissue filling each neuroforamen as one solid block per foramen:
the eight `T8-10 - additional_menginges-N.STL` bodies, four on each side, each
about 25 x 10 x 18 mm and each wrapped around one nerve root and its ganglion.
Unlike the canal compartments they are NOT hollow shells -- a ray through one
gives a single material interval, not an entry/exit pair around a cavity -- so
the "corridor" here is the block's own cross-section rather than a void inside
it. That is the physically right object anyway: a percutaneous DRG lead is
threaded THROUGH foraminal fat, displacing it, not dropped into an existing gap.

Note for anyone reading ../ansys/tissue_map.yaml: those blocks currently match
the `additional_menginges*` pattern and are therefore assigned Dura Mater. The
paper's compartment list calls this tissue "intraforaminal tissue", which is not
dura. That is a tissue-assignment question, not a geometry one, and nothing here
depends on the answer -- but it is worth settling before a production solve.

THE MEASUREMENT
---------------
At each lateral station u the block's centre is found by ALTERNATING ray casts,
seeded on the ganglion centroid and iterated to a fixed point:

    cast +Z at (x, y)  -> the interval containing z  -> new z, and the z extent
    cast +Y at (x, z)  -> the interval containing y  -> new y, and the y extent

Four rounds is more than enough; the pair converges in two. Taking the interval
that CONTAINS the current estimate rather than the first or the widest is what
keeps the walk on the ganglion's own foramen when a ray clips a neighbouring
structure.

USABLE RANGE, AND WHY IT IS NOT THE WHOLE BLOCK
-----------------------------------------------
Both ends of every block taper to nothing -- medially where it meets the canal,
laterally where the root leaves the model. In the taper the "centre" is set by
the shape of the taper rather than by any channel, and it swings by millimetres
over a couple of millimetres of u. So the usable range is the contiguous run
around the ganglion where BOTH cross-section extents are at least
MIN_EXTENT_MM = 5.0 mm, and only that run is fitted. It comes out about 20 mm
wide on every one of the eight, which is roughly two DRG contact arrays.

FIT DEGREE 4, AND IT IS NOT AS GOOD AS THE CANAL FIT
----------------------------------------------------
The dorsal canal centreline is a cubic good to 0.018 mm RMS over 98 mm. Nothing
here is that clean: the foramen is a short irregular channel between two
pedicles, and the best a quartic does is about 0.13 mm RMS in y and 0.08-0.15 mm
in z, with worst-case residuals near 0.5 mm. Degree 3 is visibly worse in z
(rms 0.3-0.46 mm, max 1.5 mm) because the z centreline rises to a peak just
outside the canal and then falls, which a cubic cannot hold; degrees 5 and 6 buy
almost nothing in y. Quartic is the knee, and the residuals are reported in the
JSON so nobody has to take that on trust. Read 0.5 mm as the honest positional
uncertainty of a generated DRG lead, against a 1.3 mm lead in a channel with
5-7 mm of room -- comfortable, but an order of magnitude looser than the canal.

VALIDATION AGAINST RADO'S OWN LEAD
-----------------------------------
RADO ships a 4-contact lead sitting 11-25 mm left of midline at z = 93-97, on a
curved insulator -- a DRG lead, on the left level-3 ganglion. It is the only
independent check available, and the corridor measured here reproduces it: over
the ganglion's own span the fitted centreline runs within about 0.6 mm in y and
0.8 mm in z of the axis of RADO's insulator, drifting to 1.3/2.8 mm only at the
extreme medial end where the block flares into the canal. See
src/freecad/build_lead_config.py --compare-rado, which measures that agreement
rather than restating it.

GANGLION NAMES
--------------
L1..L4 and R1..R4, counted ROSTRAL to CAUDAL within each side from the ganglion
cores' own z -- the same derivation body_aliases.yaml uses for its labels, and
it reproduces that file's level_map exactly. They are the model's own ganglion
numbering and NOT lumbar levels: "L3" here means "left, third from the top",
which for these bodies is a THORACIC ganglion around z = 93.

Frame (src/freecad/check_laterality.py): +X anatomical LEFT, midline x = 56.60;
+Y posterior/DORSAL; +Z rostral. Millimetres throughout.

IMPORTING THIS MODULE IS SAFE
-----------------------------
src/freecad/build_lead_config.py imports corridor_at() and the fitting helpers,
so measuring a foramen has one implementation. The guard at the bottom is what
stops that import from re-measuring all eight foramina and rewriting the JSON as
a side effect -- freecadcmd imports a script rather than exec'ing it, so the
usual __main__ guard never fires. Same story as measure_corridor.py.
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
DEFAULT_JSON = os.path.join(HERE, "foraminal_corridors.json")
DEFAULT_REPORT = os.path.join(HERE, "foramen_report.txt")

# Cross-section below which the block is tapering to its edge rather than
# offering a channel; see "USABLE RANGE" above. A DRG lead is 1.2-1.3 mm, so
# 5 mm is not a clearance rule -- it is where the measurement stops meaning
# anything.
MIN_EXTENT_MM = 5.0

# The DRG "football" cores, one per ganglion. AprilRootDGR-* is the RIGHT side
# and AprilRootDGR_right-* is the LEFT one, despite the tags: see the laterality
# note at the top of body_aliases.yaml, which check_laterality.py verifies.
DRG_CORE_GLOBS = (
    ("right", "AprilRootDGR-%d neuro_dorsal _roots_football_in_middle-1"),
    ("left", "AprilRootDGR_right-%d neuro_dorsal _roots_football_in_middle_R-1"),
)

# The eight foraminal tissue blocks. Which block belongs to which ganglion is
# NOT hardcoded -- it is found by asking which block's bounding box contains the
# ganglion core, which is the same kind of geometry-derived answer as the level
# numbering. See assign_blocks().
BLOCK_STEM = "additional_menginges-%d"
BLOCK_INDICES = tuple(range(2, 10))

# Bodies a lead must not pass through, other than the target ganglion itself:
# every nerve-root core near the corridor. Matched by filename prefix because
# tissue_map.yaml's nerve_root patterns are the authority on what is nerve and
# these are exactly its "Inside"/core members.
NEURAL_PREFIXES = ("AprilRootMainRoot", "AprilRootBranchRoot", "Assem1SmallRoot",
                   "smallRoot", "AprilRootDGR")


def stl(name):
    """Full path of one RADO STL from its stem, with or without the T8-10 prefix."""
    for candidate in ("T8-10 - %s.STL" % name, "T8-10 - %s.stl" % name,
                      "%s.stl" % name, "%s.STL" % name):
        path = os.path.join(STL_DIR, candidate)
        if os.path.exists(path):
            return path
    raise IOError("no STL for %r in %s" % (name, STL_DIR))


def polyfit(xs, ys, deg):
    """Least-squares polynomial fit, coefficients highest power first.

    Deliberately a copy of measure_corridor.polyfit rather than an import: this
    module has to be importable by build_lead_config.py without dragging in
    measure_corridor's module-level work, and the function is fifteen lines of
    normal equations with no numpy available inside freecadcmd.
    """
    n = deg + 1
    A = [[sum(x ** (i + j) for x in xs) for j in range(n)] for i in range(n)]
    b = [sum(y * x ** i for x, y in zip(xs, ys)) for i in range(n)]
    for c in range(n):
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
    return list(reversed(coef))


def polyval(coef, x):
    v = 0.0
    for c in coef:
        v = v * x + c
    return v


def area_centroid(mesh):
    """Area-weighted centroid of a mesh, as (x, y, z).

    Area weighting, not the mean of the vertices: RADO's meshes are far more
    finely tessellated where they curve, so a plain vertex mean is pulled toward
    the curved end of every body.
    """
    sx = sy = sz = sa = 0.0
    for f in mesh.Facets:
        a = f.Area
        p = f.Points
        sx += a * (p[0][0] + p[1][0] + p[2][0]) / 3.0
        sy += a * (p[0][1] + p[1][1] + p[2][1]) / 3.0
        sz += a * (p[0][2] + p[1][2] + p[2][2]) / 3.0
        sa += a
    return sx / sa, sy / sa, sz / sa


def ray_intervals(mesh, point, direction):
    """Material intervals along a ray, as [(t_in, t_out), ...] in mm from `point`.

    The signed-distance form of measure_corridor.intervals(), so it can be cast
    along +Z as well as +Y. foraminate() returns the whole infinite line, and
    sorting by the signed parameter is what pairs entries with exits correctly
    in both directions.
    """
    hits = mesh.foraminate(tuple(point), tuple(direction))
    ts = sorted((v[0] - point[0]) * direction[0]
                + (v[1] - point[1]) * direction[1]
                + (v[2] - point[2]) * direction[2] for v in hits.values())
    return [(ts[i], ts[i + 1]) for i in range(0, len(ts) - 1, 2)]


def containing(segments, want):
    """The interval containing `want`, else the widest, else None.

    Taking the containing one matters: near a ganglion a +Z ray can clip the
    neighbouring foramen's block, and "the first interval" or "the widest one"
    would then walk the centreline off to the wrong level.
    """
    for lo, hi in segments:
        if lo - 1e-9 <= want <= hi + 1e-9:
            return lo, hi
    return max(segments, key=lambda s: s[1] - s[0]) if segments else None


def corridor_at(block, x, y_seed, z_seed, rounds=4):
    """Centre and extents of the foraminal block at one x. (y, z, y_ext, z_ext).

    None when either ray misses the block entirely, which is how the caller
    learns it has walked off the end. See "THE MEASUREMENT" above for why this
    alternates rather than doing one pass.
    """
    y, z = y_seed, z_seed
    y_ext = z_ext = 0.0
    for _ in range(rounds):
        sz = containing(ray_intervals(block, (x, y, 0.0), (0.0, 0.0, 1.0)), z)
        if sz is None:
            return None
        z, z_ext = 0.5 * (sz[0] + sz[1]), sz[1] - sz[0]
        sy = containing(ray_intervals(block, (x, 0.0, z), (0.0, 1.0, 0.0)), y)
        if sy is None:
            return None
        y, y_ext = 0.5 * (sy[0] + sy[1]), sy[1] - sy[0]
    return y, z, y_ext, z_ext


def scan_corridor(block, x_sign, y_seed, z_seed, u_lo, u_hi, step):
    """Walk corridor_at() across the block. Returns [(u, y, z, y_ext, z_ext)].

    u is the lateral distance from the midline, always positive; x_sign turns it
    into a coordinate (+1 for a left-side ganglion, -1 for a right-side one), so
    one set of formulae and one sign of u describe both sides.
    """
    rows = []
    u = u_lo
    while u <= u_hi + 1e-9:
        got = corridor_at(block, MIDLINE_X + x_sign * u, y_seed, z_seed)
        if got is not None:
            rows.append((u,) + got)
        u += step
    return rows


def usable_span(rows, u_ganglion, min_extent=MIN_EXTENT_MM):
    """The contiguous run around the ganglion where the block is a real channel.

    Anchored on the ganglion rather than on the widest station, because the
    question a DRG lead asks is "how far can I go either side of THIS ganglion",
    and a block that happened to be widest somewhere else would give a run that
    does not contain the target at all.
    """
    if not rows:
        return []
    seed = min(range(len(rows)), key=lambda i: abs(rows[i][0] - u_ganglion))
    lo = seed
    while lo > 0 and rows[lo - 1][3] >= min_extent and rows[lo - 1][4] >= min_extent:
        lo -= 1
    hi = seed
    while hi < len(rows) - 1 and rows[hi + 1][3] >= min_extent \
            and rows[hi + 1][4] >= min_extent:
        hi += 1
    return rows[lo:hi + 1]


def fit_corridor(rows, degree):
    """Fit y(u) and z(u) over `rows`. Returns a dict ready for the JSON."""
    us = [r[0] for r in rows]
    out = {}
    for key, col in (("y", 1), ("z", 2)):
        coef = polyfit(us, [r[col] for r in rows], degree)
        resid = [r[col] - polyval(coef, r[0]) for r in rows]
        out["%s_poly_coef_high_to_low" % key] = coef
        out["%s_fit_rms_mm" % key] = (sum(v * v for v in resid) / len(resid)) ** 0.5
        out["%s_fit_max_abs_mm" % key] = max(abs(v) for v in resid)
    return out


def point_on(corridor, u):
    """The corridor centreline at lateral distance u, as (x, y, z).

    THE one place the stored polynomials are turned back into a point, so
    build_lead_config.py sweeping a lead and this module reporting a fit cannot
    disagree about what the corridor means.
    """
    return (MIDLINE_X + corridor["x_sign"] * u,
            polyval(corridor["y_poly_coef_high_to_low"], u),
            polyval(corridor["z_poly_coef_high_to_low"], u))


# --------------------------------------------------------------------------
# identifying which body is which -- all of it from the geometry
# --------------------------------------------------------------------------
def find_ganglia():
    """The eight DRG cores as {tag: info}, named L1..L4 / R1..R4 rostral-caudal.

    The level is the core's rank by z WITHIN ITS SIDE, which is how
    body_aliases.yaml derives the numbering in its labels; running it here from
    the geometry rather than copying that file's level_map means the two cannot
    drift apart, and it reproduces it exactly (source index 1,2,3,4 -> level
    2,4,3,1 on both sides).
    """
    found = {}
    for side, pattern in DRG_CORE_GLOBS:
        cores = []
        for i in (1, 2, 3, 4):
            path = stl(pattern % i)
            mesh = Mesh.Mesh(path)
            bb = mesh.BoundBox
            cores.append({
                "stl": os.path.basename(path),
                "source_index": i,
                "centroid": list(area_centroid(mesh)),
                "u_range": sorted([abs(bb.XMin - MIDLINE_X), abs(bb.XMax - MIDLINE_X)]),
                "z_range": [bb.ZMin, bb.ZMax],
            })
        cores.sort(key=lambda c: -c["centroid"][2])          # rostral first
        for level, core in enumerate(cores, start=1):
            core["side"] = side
            core["x_sign"] = 1 if side == "left" else -1
            found["%s%d" % (side[0].upper(), level)] = core
    return found


def assign_blocks(ganglia):
    """Attach each ganglion its foraminal block, by containment. Returns warnings.

    Asking which block's bounding box holds the ganglion core beats a hardcoded
    table for the usual reason: the answer stays right if RADO renumbers the
    files, and it fails loudly rather than silently pairing a lead with the
    wrong foramen.
    """
    blocks = []
    for i in BLOCK_INDICES:
        path = stl(BLOCK_STEM % i)
        bb = Mesh.Mesh(path).BoundBox
        blocks.append({"stl": os.path.basename(path), "index": i, "bb": bb})
    warnings = []
    for tag, core in sorted(ganglia.items()):
        cx, cy, cz = core["centroid"]
        hits = [b for b in blocks
                if b["bb"].XMin <= cx <= b["bb"].XMax
                and b["bb"].YMin <= cy <= b["bb"].YMax
                and b["bb"].ZMin <= cz <= b["bb"].ZMax]
        if len(hits) != 1:
            warnings.append("%s: %d foraminal blocks contain the ganglion core (%s)"
                            % (tag, len(hits), ", ".join(h["stl"] for h in hits) or "none"))
            core["block_stl"] = hits[0]["stl"] if hits else None
        else:
            core["block_stl"] = hits[0]["stl"]
    used = [c["block_stl"] for c in ganglia.values() if c["block_stl"]]
    if len(set(used)) != len(used):
        warnings.append("two ganglia were given the same foraminal block: %s"
                        % ", ".join(sorted(used)))
    return warnings


def neural_neighbours(corridor_bb):
    """Nerve-root and ganglion cores whose bounding box overlaps a corridor.

    Cached in the JSON so the fit check does not have to bounding-box 245 STLs
    every time it validates a lead. Bounding boxes only -- this is a shortlist
    of bodies worth testing properly, not a containment answer.
    """
    out = []
    for fn in sorted(os.listdir(STL_DIR)):
        if not fn.lower().endswith(".stl"):
            continue
        stem = fn[len("T8-10 - "):] if fn.startswith("T8-10 - ") else fn
        if not stem.startswith(NEURAL_PREFIXES):
            continue
        if "Menging" in stem or "menging" in stem or "csf" in stem:
            continue                      # the sheath shells, not the nerve itself
        bb = Mesh.Mesh(os.path.join(STL_DIR, fn)).BoundBox
        if (bb.XMin <= corridor_bb[1] and bb.XMax >= corridor_bb[0]
                and bb.YMin <= corridor_bb[3] and bb.YMax >= corridor_bb[2]
                and bb.ZMin <= corridor_bb[5] and bb.ZMax >= corridor_bb[4]):
            out.append(fn)
    return out


def main():
    argv = shlex.split(os.environ.get("MEASURE_FORAMEN_ARGS", ""))

    def opt(name, cast, default):
        return cast(argv[argv.index(name) + 1]) if name in argv else default

    step = opt("--step", float, 0.5)
    degree = opt("--degree", int, 4)
    u_lo = opt("--u-min", float, 2.0)
    u_hi = opt("--u-max", float, 34.0)

    lines = ["foraminal corridors: centre of each intraforaminal tissue block,",
             "as a function of lateral distance u from the midline (x = %.2f)."
             % MIDLINE_X, ""]

    ganglia = find_ganglia()
    for w in assign_blocks(ganglia):
        lines.append("WARNING: %s" % w)

    out = {"frame": {"midline_x": MIDLINE_X, "x_is": "anatomical left (+)",
                     "y_is": "posterior/dorsal (+)", "z_is": "rostral (+)",
                     "u_is": "lateral distance from the midline, always positive"},
           "method": "alternating +Z / +Y ray casts through the intraforaminal "
                     "block, seeded on the ganglion core centroid",
           "sample_step_mm": step, "fit_degree": degree,
           "min_extent_mm": MIN_EXTENT_MM, "ganglia": {}}

    for tag in sorted(ganglia):
        core = ganglia[tag]
        if not core["block_stl"]:
            lines.append("%s: no foraminal block -- skipped" % tag)
            continue
        block = Mesh.Mesh(os.path.join(STL_DIR, core["block_stl"]))
        _cx, cy, cz = core["centroid"]
        u_ganglion = 0.5 * (core["u_range"][0] + core["u_range"][1])

        rows = scan_corridor(block, core["x_sign"], cy, cz, u_lo, u_hi, step)
        use = usable_span(rows, u_ganglion)
        if len(use) < degree + 2:
            lines.append("%s: only %d usable stations -- cannot fit" % (tag, len(use)))
            continue
        fit = fit_corridor(use, degree)

        entry = {"side": core["side"], "x_sign": core["x_sign"],
                 "block_stl": core["block_stl"], "drg_core_stl": core["stl"],
                 "ganglion_centroid": core["centroid"],
                 "ganglion_u": core["u_range"],
                 "ganglion_u_centre": u_ganglion,
                 "usable_u": [use[0][0], use[-1][0]],
                 "n_stations": len(use),
                 "min_y_extent_mm": min(r[3] for r in use),
                 "min_z_extent_mm": min(r[4] for r in use),
                 "stations": [[round(v, 6) for v in r] for r in use]}
        entry.update(fit)

        xs = [MIDLINE_X + core["x_sign"] * u for u, _y, _z, _ye, _ze in use]
        entry["neural_stls"] = neural_neighbours(
            [min(xs) - 3, max(xs) + 3,
             min(r[1] for r in use) - 6, max(r[1] for r in use) + 6,
             min(r[2] for r in use) - 9, max(r[2] for r in use) + 9])
        out["ganglia"][tag] = entry

        lines += ["", "%s  (%s side, ganglion core at z = %.1f, %s)"
                  % (tag, core["side"], cz, core["stl"]),
                  "   foraminal block   %s" % core["block_stl"],
                  "   ganglion spans    u %.1f .. %.1f mm from the midline"
                  % (core["u_range"][0], core["u_range"][1]),
                  "   usable corridor   u %.1f .. %.1f mm (%d stations at %.2f mm)"
                  % (use[0][0], use[-1][0], len(use), step),
                  "   cross-section     at least %.2f mm dorso-ventral, %.2f mm "
                  "rostro-caudal" % (entry["min_y_extent_mm"], entry["min_z_extent_mm"]),
                  "   centreline fit    degree %d: y rms %.3f mm (max %.3f), "
                  "z rms %.3f mm (max %.3f)"
                  % (degree, fit["y_fit_rms_mm"], fit["y_fit_max_abs_mm"],
                     fit["z_fit_rms_mm"], fit["z_fit_max_abs_mm"]),
                  "   neural bodies near the corridor: %d"
                  % len(entry["neural_stls"])]

    json.dump(out, open(DEFAULT_JSON, "w"), indent=2)
    lines += ["", "wrote %s (%d corridors)"
              % (os.path.relpath(DEFAULT_JSON, os.path.dirname(HERE)),
                 len(out["ganglia"]))]
    open(DEFAULT_REPORT, "w").write("\n".join(lines) + "\n")
    sys.stdout.write("\n".join(lines) + "\n")


def run_as_freecadcmd_script():
    """True when a FreeCAD interpreter was handed THIS file to run.

    Same guard as measure_corridor.py, and for the same reason: without it,
    `from measure_foramen import corridor_at` in build_lead_config.py would
    re-measure all eight foramina and rewrite foraminal_corridors.json as an
    import side effect.
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
