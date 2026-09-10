#!/usr/bin/env python3
"""Build and VALIDATE any named SCS lead configuration from lead_configs.yaml.

WHAT A CLINICIAN ASKS FOR, AND WHAT THIS TURNS IT INTO
------------------------------------------------------
The clinical description of a lead is "8 contacts, 3 mm long, 1 mm apart,
dorsal, at T10, 1 mm left of midline". None of that is geometry. This module is
the layer that turns it into geometry and then, more importantly, tells you
whether the answer is anatomically possible in THIS model:

    lead_configs.yaml  ->  resolve()   ->  a flat set of numbers
                       ->  validate()  ->  PASS/WARN/FAIL with the measurements
                       ->  export()    ->  generated_leads/<name>/*.stl
                       ->  preview()   ->  disposable bodies in the open document

src/freecad/lead_designer.FCMacro is the same thing with spinboxes on it.

WHY THERE IS NO BOOLEAN SUBTRACTION HERE
----------------------------------------
The intuitive mental model is that a lead has to be CARVED out of the tissue
around it -- Mesh.difference() the epidural fat against the lead so the two do
not occupy the same space. That is not what this model needs, for two reasons,
and the fit check below is what replaces it.

1. RADO's compartments are HOLLOW SHELLS that tile space. Each body already has
   the inner ones carved out of it, and the epidural fat is a real void that the
   lead sits INSIDE -- the same void a percutaneous lead is threaded into in
   theatre. Ray casting the epidural mesh (see measure_corridor.py) measures
   that void directly: 2.27-2.35 mm dorsally, 1.67-1.74 mm ventrally at the
   midline. A 1.3 mm lead fits either side with clearance to spare. There is
   nothing to subtract: the lead displaces fat that the mesher will simply not
   fill, exactly as the FEM wants.

2. Mesh.difference() is a COMPLETE NO-OP in this FreeCAD build -- verified. It
   returns without error and changes nothing. Anything built on it would look
   like it worked and quietly produce an unmodified model.

So the correct check is not "did the subtraction succeed" but "is carving
unnecessary", i.e. does the lead lie wholly within the fat void with clearance.
That is what validate() answers, in three measurements:

    fit         the channel thickness measured AT the requested side, lateral
                offset and z span, against the lead diameter
    containment every sampled vertex of the built lead tested for being inside
                the epidural fat and outside the dura -- the same check that
                verified the existing dorsal/ventral leads at 100% / 0%
    coverage    whether the lead runs off the end of the measured centreline

A configuration that fails is REPORTED, in clinical terms, and not exported. It
is not carved into the anatomy to make it fit. If a lead does not fit the
epidural space, the honest answer is that it does not fit the epidural space.

WHY LATERAL OFFSETS RE-MEASURE THE CENTRELINE
---------------------------------------------
epidural_corridor.json holds the centreline y(z) measured at the MIDLINE. The
channel centre moves as you go lateral -- about 1.1 mm of y at 4 mm off the
midline -- so a lead swept along the midline curve at x = midline + 4 sits about
half its clearance off centre before it has done anything else. Any configuration
with a non-zero x_offset therefore gets its own centreline, ray cast and fitted
at ITS x over ITS z span. At x_offset = 0 this reproduces the stored polynomial;
the report prints the agreement so the two can never silently diverge.

The sweep itself is NOT reimplemented here. make_scs_lead.segment_plan() and
make_scs_lead.sweep_segments() are imported, so "sweep a lead along the measured
canal centreline" has exactly one implementation in the repo; this module only
decides which curve to hand it. Likewise the ray cast and the polynomial fit come
from measure_corridor.

NO .FCStd PER CONFIGURATION
---------------------------
The document is 25 MB. Ten configurations must not mean ten documents, and every
save of one is another full copy in git history. So the durable artefacts are
this YAML entry plus generated_leads/<name>/, both small and both text-ish, and
the bodies in a FreeCAD document are a DISPOSABLE PREVIEW:

    preview(doc, params)   adds the lead to the open document inside a group
                           named SCS_LeadPreview, replacing any previous preview
    clear_preview(doc)     deletes that group and everything in it

Explore ten configurations, keep the one you want, export its STLs, and close
the document WITHOUT SAVING. The two committed variants
(NBF_RADO-SCS_dorsal.FCStd / _ventral.FCStd) exist because they are the
dorsal-vs-ventral study's two arms and are handed to Ansys; they were built by
make_lead_variants.py, and they are the only two that should ever exist.

NEVER SAVE A DOCUMENT FROM HERE
-------------------------------
This module contains no doc.save() and must not grow one. A document saved
without a running GUI is written with no GuiDocument.xml and no ShapeAppearance
blobs at all -- every tissue colour gone, the document still opening and looking
like a fresh grey import. make_tissue_groups.py documents the measurement. The
preview path additionally REFUSES to touch a document that has a FileName when
FreeCAD.GuiUp is false, so a headless run cannot leave a real model dirty and
tempt someone into saving it.

USAGE
-----
The listing and resolving half needs no FreeCAD at all:

    python3 src/freecad/build_lead_config.py --list
    python3 src/freecad/build_lead_config.py --name khadka_clinical_T10 --show

The measuring and building half does:

    BUILD_LEAD_CONFIG_ARGS="--name khadka_clinical_T10 --validate" \\
        freecadcmd src/freecad/build_lead_config.py
    BUILD_LEAD_CONFIG_ARGS="--all --validate" \\
        freecadcmd src/freecad/build_lead_config.py
    BUILD_LEAD_CONFIG_ARGS="--name dorsal_T11_left1 --export" \\
        freecadcmd src/freecad/build_lead_config.py

Any parameter can be overridden on top of a named entry, which is how the macro
drives it for a custom lead:

    BUILD_LEAD_CONFIG_ARGS="--name dorsal_z110_8c --level T10 --x-offset -2 --validate"

Options come from BUILD_LEAD_CONFIG_ARGS under any FreeCAD interpreter, because
FreeCAD owns the command line; output goes to a report file because FreeCAD
swallows print(). Same convention as every other script here.
"""
import argparse
import json
import os
import re
import shlex
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
DEFAULT_CONFIGS = os.path.join(HERE, "lead_configs.yaml")
DEFAULT_CORRIDOR = os.path.join(HERE, "epidural_corridor.json")
DEFAULT_STL = os.path.join(REPO, "STL_files")
DEFAULT_OUT = os.path.join(HERE, "generated_leads")
DEFAULT_REPORT = os.path.join(HERE, "build_lead_config_report.txt")

MIDLINE_X = 56.60
EPIDURAL_STL = "T8-10 - neuro_EpiduralSpace-1.STL"
DURA_STL = "T8-10 - neuro_Meninges-1.STL"

# The name of the group every previewed body lives in. Matched on Name, not
# Label, so retitling it in the tree does not orphan the preview.
PREVIEW_GROUP = "SCS_LeadPreview"

# Clearance below which a fit is reported as tight rather than comfortable. Half
# of it is the amount the lead can wander off the channel centre before it
# touches a wall, and the centreline fit itself is good to about 0.02 mm, so
# 0.20 mm is roughly "ten times the modelling error" rather than a clinical rule.
TIGHT_CLEARANCE_MM = 0.20

# Every parameter a configuration may carry, and how to read it. "level" and
# "note" are strings; everything else is a number. Anything not in here is
# rejected by resolve() rather than silently ignored -- a typo'd "diamter: 2.0"
# that quietly kept the 1.3 mm default would be the worst kind of bug here.
PARAM_TYPES = {
    "side": str, "contacts": int, "contact_length": float, "gap": float,
    "diameter": float, "x_offset": float, "z_center": float, "tail": float,
    "level": str, "note": str,
}


def is_freecad_interpreter():
    """True under freecadcmd or the GUI console; see apply_colors.py."""
    return os.path.basename(sys.argv[0]).lower().startswith("freecad")


def script_args():
    """Arguments meant for this script, under either interpreter."""
    if is_freecad_interpreter():
        return shlex.split(os.environ.get("BUILD_LEAD_CONFIG_ARGS", ""))
    return sys.argv[1:]


# --------------------------------------------------------------------------
# lead_configs.yaml
# --------------------------------------------------------------------------
def parse_config_file(path):
    """Read lead_configs.yaml. Returns {"defaults": {}, "levels": {}, "configs": {}}.

    A hand-rolled parser, for the reason requirements.txt gives: PyYAML is
    deliberately not a dependency of this directory, so the scripts run on a bare
    interpreter with nothing installed. It accepts exactly the subset
    lead_configs.yaml uses and nothing else -- a top-level mapping whose values
    are mappings of scalars, at most two levels deep, with # comments. No lists,
    no flow style, no anchors. Anything it cannot parse raises rather than being
    skipped, because a silently dropped parameter is a wrong lead.
    """
    top, section, entry = {}, None, None
    for n, raw in enumerate(open(path), 1):
        line = raw.rstrip("\n")
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        indent = len(line) - len(line.lstrip())
        m = re.match(r"^\s*([A-Za-z_][A-Za-z0-9_.\-]*)\s*:\s*(.*)$", line)
        if not m:
            raise ValueError("%s:%d: cannot parse %r" % (path, n, line))
        key, rest = m.group(1), strip_comment(m.group(2))

        if indent == 0:
            if rest:
                raise ValueError("%s:%d: top-level %r must be a mapping" % (path, n, key))
            section, entry = {}, None
            top[key] = section
        elif indent == 2:
            if section is None:
                raise ValueError("%s:%d: %r is outside any section" % (path, n, key))
            if rest == "":
                entry = {}
                section[key] = entry
            else:
                section[key] = scalar(rest)
                entry = None
        elif indent == 4:
            if entry is None:
                raise ValueError("%s:%d: %r is outside any entry" % (path, n, key))
            entry[key] = scalar(rest)
        else:
            raise ValueError("%s:%d: unexpected indent %d" % (path, n, indent))

    for required in ("defaults", "levels", "configs"):
        if required not in top:
            raise ValueError("%s: no %r section" % (path, required))
    return top


def strip_comment(text):
    """Drop a trailing # comment, unless it is inside quotes."""
    out, quote = [], None
    for i, ch in enumerate(text):
        if quote:
            out.append(ch)
            if ch == quote:
                quote = None
        elif ch in "\"'":
            quote = ch
            out.append(ch)
        elif ch == "#":
            break
        else:
            out.append(ch)
    return "".join(out).strip()


def scalar(text):
    """int, float or (unquoted) string."""
    if len(text) >= 2 and text[0] == text[-1] and text[0] in "\"'":
        return text[1:-1]
    try:
        return int(text)
    except ValueError:
        pass
    try:
        return float(text)
    except ValueError:
        pass
    return text


def resolve(name, cfg, overrides=None):
    """Flatten one named configuration into the numbers the geometry needs.

    Layering, later wins: defaults -> the named entry -> overrides (the CLI
    flags, or the macro's spinboxes). `level: T10` is turned into the z of that
    vertebral body here, so nothing downstream has to know what a level is; an
    explicit z_center in a LATER layer beats an earlier level, and vice versa,
    which is what lets the macro offer both controls without them fighting.

    Returns a plain dict, plus "name", "x" (absolute) and "level" for reporting.
    """
    if name is not None and name not in cfg["configs"]:
        raise ValueError("no configuration named %r (have: %s)"
                         % (name, ", ".join(sorted(cfg["configs"]))))

    params = dict(cfg["defaults"])
    order = ["defaults"]
    if name is not None:
        params.update(cfg["configs"][name])
        order.append(name)
    if overrides:
        params.update({k: v for k, v in overrides.items() if v is not None})
        order.append("overrides")

    unknown = [k for k in params if k not in PARAM_TYPES]
    if unknown:
        raise ValueError("unknown parameter(s) %s in configuration %r"
                         % (", ".join(sorted(unknown)), name))
    for k, v in list(params.items()):
        if isinstance(v, str) and PARAM_TYPES[k] is not str:
            params[k] = PARAM_TYPES[k](v)

    level = params.get("level")
    if level:
        if level not in cfg["levels"]:
            raise ValueError("no vertebral level named %r (have: %s)"
                             % (level, ", ".join(cfg["levels"])))
        # A level named in the same or a later layer than an explicit z_center
        # wins; that is what "later wins" means when both are present, and the
        # only way to tell is which layer each came from.
        z_from_level = float(cfg["levels"][level])
        if _layer_of("z_center", cfg, name, overrides, order) <= \
           _layer_of("level", cfg, name, overrides, order):
            params["z_center"] = z_from_level

    params["name"] = name or "custom"
    params["x"] = MIDLINE_X + float(params["x_offset"])
    check_sane(params)
    return params


def _layer_of(key, cfg, name, overrides, order):
    """Index of the last layer that set `key`; -1 if none did. See resolve()."""
    layers = [cfg["defaults"]]
    if name is not None:
        layers.append(cfg["configs"][name])
    if overrides:
        layers.append({k: v for k, v in overrides.items() if v is not None})
    last = -1
    for i, layer in enumerate(layers):
        if key in layer:
            last = i
    return last


def check_sane(params):
    """Reject nonsense before any geometry is attempted, in the user's terms."""
    problems = []
    if params["side"] not in ("dorsal", "ventral"):
        problems.append("side must be dorsal or ventral, not %r" % params["side"])
    if params["contacts"] < 1:
        problems.append("a lead needs at least one contact")
    for key, what in (("contact_length", "contact length"), ("diameter", "diameter")):
        if params[key] <= 0:
            problems.append("%s must be greater than zero" % what)
    if params["gap"] < 0:
        problems.append("the gap between contacts cannot be negative")
    if params["tail"] < 0:
        problems.append("the tail cannot be negative")
    if problems:
        raise ValueError("; ".join(problems))


def array_length(params):
    """Length of the contact array, mm (first contact's start to last one's end)."""
    return (params["contacts"] * params["contact_length"]
            + (params["contacts"] - 1) * params["gap"])


def describe(params):
    """The configuration as a clinician would say it out loud."""
    side_word = {"dorsal": "dorsal", "ventral": "ventral"}[params["side"]]
    dx = float(params["x_offset"])
    if abs(dx) < 1e-9:
        where = "on the midline"
    else:
        where = "%.1f mm %s of midline" % (abs(dx), "left" if dx > 0 else "right")
    at = params.get("level") or "z = %.1f mm" % params["z_center"]
    return ("%d contacts, %.1f mm long, %.1f mm apart, %.2f mm diameter, "
            "in the %s epidural space at %s, %s"
            % (params["contacts"], params["contact_length"], params["gap"],
               params["diameter"], side_word, at, where))


# --------------------------------------------------------------------------
# measurement and geometry -- everything below needs FreeCAD
# --------------------------------------------------------------------------
def _geometry_modules():
    """Import the FreeCAD-only modules. Deferred so --list/--show run anywhere.

    make_scs_lead and measure_corridor both import FreeCAD at module scope, and
    both call main() at module scope UNLESS handed to freecadcmd as the script --
    which is why they each grew a run_as_freecadcmd_script() guard when this
    module started importing them. Without it, importing measure_corridor here
    would silently re-measure the corridor and rewrite epidural_corridor.json.
    """
    if HERE not in sys.path:
        sys.path.insert(0, HERE)
    import FreeCAD
    import Mesh
    import make_scs_lead
    import measure_corridor
    return FreeCAD, Mesh, make_scs_lead, measure_corridor


def load_meshes(stl_dir=DEFAULT_STL):
    """The two meshes every check needs: the epidural fat and the dura."""
    _fc, Mesh, _msl, _mc = _geometry_modules()
    return {"epidural": Mesh.Mesh(os.path.join(stl_dir, EPIDURAL_STL)),
            "dura": Mesh.Mesh(os.path.join(stl_dir, DURA_STL))}


def channel_at(epi, x, z, side):
    """The epidural channel on one side at (x, z): (y_lo, y_hi, thickness).

    None when the ray does not give exactly two material intervals. Two is the
    signature of a clean slice -- ventral fat, then the thecal sac, then dorsal
    fat -- and anything else means the ray went through a foramen, a root sleeve
    or off the end of the sac, where "the thickness of the channel" is not a
    well-defined quantity. Reported rather than averaged over.
    """
    _fc, _Mesh, _msl, mc = _geometry_modules()
    segs = mc.intervals(epi, x, z)
    if len(segs) != 2:
        return None
    lo, hi, thick = segs[1] if side == "dorsal" else segs[0]
    return lo, hi, thick


def measure_centreline(epi, x, side, z_lo, z_hi, step=1.0, degree=3):
    """Ray cast the channel at this x over this z span and fit its centre y(z).

    Returns (coef, stations, skipped) where stations is [(z, y_centre, thickness)]
    and skipped is the list of z where the slice was not clean. The fit is the
    same cubic-by-normal-equations measure_corridor.py uses; the sampling is
    finer (1 mm rather than 2 mm) because it only has to cover one lead rather
    than the whole model.

    Raises when there are too few clean slices to fit -- that is a real answer
    about the anatomy at that x, not a numerical inconvenience, and the caller
    turns it into a failed validation.
    """
    _fc, _Mesh, _msl, mc = _geometry_modules()
    bb = epi.BoundBox
    z = max(z_lo, bb.ZMin + 0.5)
    stop = min(z_hi, bb.ZMax - 0.5)
    stations, skipped = [], []
    while z <= stop + 1e-9:
        seg = channel_at(epi, x, z, side)
        if seg is None:
            skipped.append(z)
        else:
            stations.append((z, 0.5 * (seg[0] + seg[1]), seg[2]))
        z += step
    if len(stations) < degree + 2:
        raise ValueError(
            "only %d clean slices of the %s channel at x = %.2f over z %.1f..%.1f -- "
            "not enough to follow it" % (len(stations), side, x, z_lo, z_hi))
    coef = mc.polyfit([s[0] for s in stations], [s[1] for s in stations], degree)
    return coef, stations, skipped


def centreline_for(params, epi, corridor, span):
    """Choose the curve to sweep along, and say why. Returns (y_of_z, info).

    At the midline the stored polynomial from epidural_corridor.json is used --
    it is measured over the whole 98 mm at 2 mm and is the curve the existing
    generated leads were built on, so reusing it keeps them bit-for-bit
    reproducible. Off the midline it is the wrong curve (see the module
    docstring), so the channel is re-measured at the requested x.

    info carries the agreement between the two at the midline, which is a live
    check that the local measurement and the stored fit have not diverged.
    """
    _fc, _Mesh, _msl, mc = _geometry_modules()
    side = params["side"]
    stored = corridor["sides"][side]["centre_poly_coef_high_to_low"]
    z_lo, z_hi = span

    local_coef, stations, skipped = measure_centreline(
        epi, params["x"], side, z_lo, z_hi)
    diff = [abs(mc.polyval(stored, z) - mc.polyval(local_coef, z))
            for z, _y, _t in stations]
    info = {"stations": stations, "skipped": skipped,
            "max_disagreement_mm": max(diff) if diff else 0.0}

    if abs(float(params["x_offset"])) < 1e-9:
        info["source"] = "stored midline centreline (epidural_corridor.json)"
        coef = stored
    else:
        info["source"] = ("centreline re-measured at x = %.2f (%d slices)"
                          % (params["x"], len(stations)))
        coef = local_coef
    info["coef"] = coef
    return (lambda z: mc.polyval(coef, z)), info


def build_shapes(params, y_of_z):
    """The lead itself: [(kind, index, shape)], plus the segment plan.

    One line of geometry, and it is not written here -- see the module docstring
    on why make_scs_lead owns the sweep.
    """
    _fc, _Mesh, msl, _mc = _geometry_modules()
    segments, z0, total = msl.segment_plan(
        params["contacts"], params["contact_length"], params["gap"],
        params["tail"], params["z_center"])
    shapes = msl.sweep_segments(segments, params["x"], y_of_z,
                                0.5 * params["diameter"])
    return shapes, segments, z0, total


def sample_points(shapes, tolerance=0.02, max_samples=6000):
    """Vertices of the built lead, for the containment check.

    Tessellating each segment puts vertices on the two end circles of every
    cylinder -- which is also all an exported STL of a cylinder has, so this is
    the same point set the check that verified the existing leads ran on. The
    straight sides mean nothing samples the middle of a segment, and that is
    fine here: a segment is at most `tail` mm long, and the centreline's
    curvature (y'' ~= 0.0035 /mm) puts the mid-segment excursion at
    y''*L^2/8 ~= 0.016 mm for a 6 mm tail -- two orders below the clearance
    being measured.

    Subsampled by a fixed stride rather than randomly, so the same configuration
    always reports the same numbers.
    """
    pts = []
    for _kind, _idx, shape in shapes:
        verts = shape.tessellate(tolerance)[0]
        pts += [(v.x, v.y, v.z) for v in verts]
    if len(pts) > max_samples:
        stride = len(pts) // max_samples + 1
        pts = pts[::stride]
    return pts


def inside(mesh, point, direction=(0.0, 1.0, 0.0)):
    """Is `point` inside the material of this mesh? Crossing parity along a ray.

    foraminate() returns the intersections of the whole INFINITE line, not the
    forward ray, so the crossings behind the point have to be dropped by hand --
    counting all of them would make every point look outside. Odd number of
    crossings ahead => inside the material. This is the test measure_corridor.py
    documents: RADO's compartments are hollow shells, so "inside the epidural
    body" means inside its fat, not inside its central cavity.
    """
    hits = mesh.foraminate(point, direction)
    n = 0
    for v in hits.values():
        t = ((v[0] - point[0]) * direction[0] + (v[1] - point[1]) * direction[1]
             + (v[2] - point[2]) * direction[2])
        if t > 1e-9:
            n += 1
    return n % 2 == 1


# --------------------------------------------------------------------------
# validation
# --------------------------------------------------------------------------
def validate(params, meshes=None, corridor=None, max_samples=6000):
    """Can this lead exist in this anatomy? Returns a result dict.

    result["checks"] is [(level, title, message)] with level in PASS/WARN/FAIL,
    written to be read by someone who has never opened FreeCAD. result["ok"] is
    False as soon as anything FAILs; a WARN is something to know about, not a
    reason to refuse -- the commonest one is a lead whose passive tail runs past
    the end of the measured corridor while every contact is comfortably inside.

    result["numbers"] carries the measurements themselves so the caller can
    print or plot them, and result["shapes"] the built geometry, so a caller
    that is about to preview or export does not rebuild it.
    """
    if meshes is None:
        meshes = load_meshes()
    if corridor is None:
        corridor = json.load(open(DEFAULT_CORRIDOR))

    epi, dura = meshes["epidural"], meshes["dura"]
    side = params["side"]
    checks, numbers = [], {}

    length = array_length(params) + 2 * params["tail"]
    z0 = params["z_center"] - 0.5 * length
    z1 = z0 + length
    c0 = params["z_center"] - 0.5 * array_length(params)
    c1 = c0 + array_length(params)
    numbers.update(lead_z=[z0, z1], contacts_z=[c0, c1], length_mm=length,
                   array_mm=array_length(params), x=params["x"])

    # -- 1. does the centreline reach this far? --------------------------
    usable = corridor["sides"][side]["usable_z"]
    numbers["usable_z"] = usable
    if c0 < usable[0] or c1 > usable[1]:
        checks.append(("FAIL", "Reaches past the model",
                       "The contacts would sit at z %.1f-%.1f mm, but this model's "
                       "%s epidural space is only measured over z %.1f-%.1f mm. "
                       "Part of the active lead would be outside the anatomy. Move "
                       "the lead %.1f mm %s, or use fewer contacts."
                       % (c0, c1, side, usable[0], usable[1],
                          max(usable[0] - c0, c1 - usable[1]),
                          "rostral" if c0 < usable[0] else "caudal")))
    elif z0 < usable[0] or z1 > usable[1]:
        checks.append(("WARN", "Tail runs off the measured corridor",
                       "Every contact is inside the measured anatomy (z %.1f-%.1f mm), "
                       "but the passive tail reaches z %.1f-%.1f mm, past the "
                       "measured range z %.1f-%.1f mm, so the last %.1f mm of "
                       "insulator follows an extrapolated curve. Harmless for the "
                       "field, worth knowing before meshing."
                       % (c0, c1, z0, z1, usable[0], usable[1],
                          max(usable[0] - z0, z1 - usable[1]))))
    else:
        checks.append(("PASS", "Inside the modelled levels",
                       "The lead spans z %.1f-%.1f mm, within the measured %s "
                       "corridor z %.1f-%.1f mm." % (z0, z1, side, usable[0], usable[1])))

    # -- 2. build it, following the right curve ---------------------------
    try:
        y_of_z, cinfo = centreline_for(params, epi, corridor, (z0, z1))
    except ValueError as exc:
        checks.append(("FAIL", "No usable channel here",
                       "The %s epidural space could not be followed at this "
                       "position: %s. This usually means the lead is off the end "
                       "of the dural sac or beside a nerve-root sleeve rather "
                       "than in the canal." % (side, exc)))
        return {"ok": False, "checks": checks, "numbers": numbers,
                "shapes": None, "params": params}

    numbers["centreline_source"] = cinfo["source"]
    numbers["centreline_vs_stored_max_mm"] = cinfo["max_disagreement_mm"]
    numbers["skipped_slices"] = len(cinfo["skipped"])
    shapes, segments, _z0, _total = build_shapes(params, y_of_z)

    if cinfo["skipped"]:
        checks.append(("WARN", "Some slices could not be measured",
                       "%d of %d sampled cross-sections along the lead did not give "
                       "a clean channel (at z %s). Those are usually foramen levels "
                       "where a nerve root leaves; the centreline is interpolated "
                       "through them."
                       % (len(cinfo["skipped"]),
                          len(cinfo["skipped"]) + len(cinfo["stations"]),
                          ", ".join("%.0f" % z for z in cinfo["skipped"][:8]))))

    # -- 3. does the lead fit the channel where it actually lies? ---------
    thicknesses = [t for _z, _y, t in cinfo["stations"]]
    t_min, t_med = min(thicknesses), sorted(thicknesses)[len(thicknesses) // 2]
    z_at_min = min(cinfo["stations"], key=lambda s: s[2])[0]
    clearance = t_min - params["diameter"]
    numbers.update(channel_min_mm=t_min, channel_median_mm=t_med,
                   channel_min_at_z=z_at_min, clearance_mm=clearance,
                   widest_lead_mm=t_min)
    where = ("on the midline" if abs(float(params["x_offset"])) < 1e-9
             else "%.1f mm %s of midline" % (abs(float(params["x_offset"])),
                                             "left" if float(params["x_offset"]) > 0
                                             else "right"))
    if clearance < 0:
        checks.append(("FAIL", "Too thick for this space",
                       "The lead is %.2f mm across, but the %s epidural fat %s is "
                       "only %.2f mm thick at its narrowest (z = %.0f mm). It would "
                       "press through the dura. The widest lead that fits here is "
                       "%.2f mm; the %s space is the roomier one at %.2f mm."
                       % (params["diameter"], side, where, t_min, z_at_min, t_min,
                          "dorsal" if side == "ventral" else "ventral",
                          corridor["sides"]["dorsal" if side == "ventral"
                                            else "ventral"]["thickness_min_mm"])))
    elif clearance < TIGHT_CLEARANCE_MM:
        checks.append(("WARN", "A tight fit",
                       "The lead is %.2f mm across and the %s fat %s is %.2f mm at "
                       "its narrowest (z = %.0f mm) -- only %.2f mm to spare, i.e. "
                       "%.2f mm each side. It fits, but it is against both walls."
                       % (params["diameter"], side, where, t_min, z_at_min,
                          clearance, 0.5 * clearance)))
    else:
        checks.append(("PASS", "Fits the epidural space",
                       "The %s fat %s is %.2f mm thick at its narrowest along this "
                       "lead (median %.2f mm) against a %.2f mm lead: %.2f mm of "
                       "clearance, %.2f mm each side."
                       % (side, where, t_min, t_med, params["diameter"],
                          clearance, 0.5 * clearance)))

    # -- 4. is the built lead actually in the fat, and out of the dura? ---
    pts = sample_points(shapes, max_samples=max_samples)
    in_fat = [p for p in pts if inside(epi, p)]
    in_dura = [p for p in pts if inside(dura, p)]
    pct_fat = 100.0 * len(in_fat) / len(pts)
    pct_dura = 100.0 * len(in_dura) / len(pts)
    numbers.update(sampled_vertices=len(pts), pct_in_fat=pct_fat, pct_in_dura=pct_dura)

    if pct_dura > 0:
        stray = [p for p in pts if inside(dura, p)]
        checks.append(("FAIL", "Overlaps the dura",
                       "%.1f%% of %d sampled points on the lead surface are inside "
                       "the dura, between z %.1f and %.1f mm. The lead is not in the "
                       "epidural space there -- it is in the thecal sac."
                       % (pct_dura, len(pts), min(p[2] for p in stray),
                          max(p[2] for p in stray))))
    elif pct_fat < 100.0:
        stray = [p for p in pts if not inside(epi, p)]
        checks.append(("FAIL", "Leaves the epidural fat",
                       "%.1f%% of %d sampled points on the lead surface are inside "
                       "the epidural fat; the other %.1f%% are outside it, between "
                       "z %.1f and %.1f mm. The lead would be sitting partly in "
                       "bone or in the foramen."
                       % (pct_fat, len(pts), 100.0 - pct_fat,
                          min(p[2] for p in stray), max(p[2] for p in stray))))
    else:
        checks.append(("PASS", "Entirely within the epidural fat",
                       "All %d sampled points on the lead surface are inside the "
                       "epidural fat and none are in the dura." % len(pts)))

    ok = not any(level == "FAIL" for level, _t, _m in checks)
    return {"ok": ok, "checks": checks, "numbers": numbers,
            "shapes": shapes, "segments": segments, "params": params}


def report_lines(params, result):
    """The validation as text, for the report file and the macro's panel."""
    n = result["numbers"]
    lines = ["%s -- %s" % (params["name"], "PASS" if result["ok"] else "FAIL"),
             describe(params), ""]
    if params.get("note"):
        lines += [params["note"], ""]
    for level, title, message in result["checks"]:
        lines.append("[%-4s] %s" % (level, title))
        for chunk in wrap(message, 74):
            lines.append("        " + chunk)
    lines.append("")
    lines.append("measurements")
    lines.append("   lead        z %.1f .. %.1f mm (%.1f mm), contacts z %.1f .. %.1f"
                 % (n["lead_z"][0], n["lead_z"][1], n["length_mm"],
                    n["contacts_z"][0], n["contacts_z"][1]))
    lines.append("   position    x = %.2f mm (midline %.2f%+.2f)"
                 % (n["x"], MIDLINE_X, float(params["x_offset"])))
    if "channel_min_mm" in n:
        lines.append("   channel     min %.2f mm at z %.1f, median %.2f mm; lead %.2f mm; "
                     "clearance %+.2f mm"
                     % (n["channel_min_mm"], n["channel_min_at_z"], n["channel_median_mm"],
                        params["diameter"], n["clearance_mm"]))
        lines.append("   containment %.1f%% of %d sampled vertices in epidural fat, "
                     "%.1f%% in dura"
                     % (n["pct_in_fat"], n["sampled_vertices"], n["pct_in_dura"]))
        lines.append("   centreline  %s" % n["centreline_source"])
        lines.append("   local fit vs stored midline fit: max %.3f mm apart"
                     % n["centreline_vs_stored_max_mm"])
    return lines


def wrap(text, width):
    """Minimal greedy wrap; textwrap would do, but this keeps the output stable."""
    words, line, out = text.split(), "", []
    for w in words:
        if line and len(line) + 1 + len(w) > width:
            out.append(line)
            line = w
        else:
            line = (line + " " + w) if line else w
    if line:
        out.append(line)
    return out


# --------------------------------------------------------------------------
# outputs: STLs on disk, and a disposable preview in a live document
# --------------------------------------------------------------------------
def export_stls(params, shapes, out_root=DEFAULT_OUT, tag=None):
    """Write one STL per contact plus one fused insulator. Returns (dir, [files]).

    Filenames are RADO's own -- "SCS Lead Electrode 1.stl", "SCS Lead
    Insulator.stl" -- deliberately: apply_colors.py and ../ansys/tissue_map.yaml
    match the object Name those sanitize to, so renaming them would cost the
    contacts their silver and the insulator its material assignment. See
    make_lead_variants.py for what FreeCAD then does to the Names on import.

    Nothing here opens or saves a project document; the scratch document exists
    only because Mesh.export() takes document objects.
    """
    FreeCAD, Mesh, msl, _mc = _geometry_modules()
    outdir = os.path.join(out_root, tag or params["name"])
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    doc = FreeCAD.newDocument("scs_lead_export")
    written = []
    try:
        for kind, idx, shape in shapes:
            if kind != "contact":
                continue
            obj = doc.addObject("Part::Feature", "contact%d" % idx)
            obj.Shape = shape
            path = os.path.join(outdir, "SCS Lead Electrode %d.stl" % idx)
            Mesh.export([obj], path)
            written.append(os.path.basename(path))
        ins = doc.addObject("Part::Feature", "insulator")
        ins.Shape = msl.fuse_insulator(shapes)
        path = os.path.join(outdir, "SCS Lead Insulator.stl")
        Mesh.export([ins], path)
        written.append(os.path.basename(path))
    finally:
        FreeCAD.closeDocument(doc.Name)
    return outdir, written


def clear_preview(doc):
    """Delete the previous preview from `doc`. Returns how many objects went.

    Found by the group's Name, which FreeCAD keeps stable, rather than its Label,
    which someone may have retitled. Children are removed explicitly because
    removing a group object on its own leaves them behind at the top of the tree.
    """
    group = doc.getObject(PREVIEW_GROUP)
    if group is None:
        return 0
    names = [o.Name for o in group.Group]
    for name in names:
        doc.removeObject(name)
    doc.removeObject(group.Name)
    doc.recompute()
    return len(names) + 1


def preview(doc, params, shapes, colours=None):
    """Put the lead into a live document as a disposable preview.

    Replaces any previous preview, so exploring configurations does not silt the
    tree up. The bodies are Part::Feature -- parametric shapes, a few kB -- not
    imported meshes, so a preview costs nothing and leaves nothing on disk. Which
    is the whole point: the configuration is the artefact, the document is not.

    THIS FUNCTION DOES NOT SAVE, and must not learn to. It refuses outright on a
    document that has a FileName when there is no GUI, because that is the shape
    of the accident that destroys the model: a headless run leaves the real
    document dirty, someone saves it, and every ShapeAppearance in it is gone.
    See make_tissue_groups.py.
    """
    FreeCAD, _Mesh, msl, _mc = _geometry_modules()
    if not FreeCAD.GuiUp and getattr(doc, "FileName", ""):
        raise RuntimeError(
            "refusing to preview into %s with no GUI running: a headless save "
            "would strip the document's colours. Preview from a running FreeCAD, "
            "or export STLs instead." % os.path.basename(doc.FileName))

    clear_preview(doc)
    group = doc.addObject("App::DocumentObjectGroup", PREVIEW_GROUP)
    group.Label = "PREVIEW: %s" % params["name"]

    made = []
    for kind, idx, shape in shapes:
        if kind != "contact":
            continue
        obj = doc.addObject("Part::Feature", "SCS_Preview_Contact%d" % idx)
        obj.Shape = shape
        obj.Label = "preview contact %d" % idx
        made.append(("electrode_contact", obj))
    ins = doc.addObject("Part::Feature", "SCS_Preview_Insulator")
    ins.Shape = msl.fuse_insulator(shapes)
    ins.Label = "preview insulator"
    made.append(("lead_insulation", ins))

    group.Group = [obj for _t, obj in made]
    if colours:
        for tissue, obj in made:
            rgb = colours.get(tissue)
            if rgb and obj.ViewObject is not None:
                paint(obj.ViewObject, rgb)
    doc.recompute()
    return group, [obj for _t, obj in made]


def paint(vo, rgb):
    """Colour one preview body. Cosmetic only -- nothing here is ever saved.

    Writes ShapeAppearance where the build has it and falls back to ShapeColor,
    the same two-step apply_colors.py uses; a preview that came up grey would be
    hard to read against the tissue colours around it.
    """
    colour = (rgb[0] / 255.0, rgb[1] / 255.0, rgb[2] / 255.0, 0.0)
    if hasattr(vo, "ShapeAppearance"):
        try:
            mat = vo.ShapeAppearance[0]
            mat.DiffuseColor = colour
            vo.ShapeAppearance = [mat]
            return
        except Exception:
            pass
    if hasattr(vo, "ShapeColor"):
        vo.ShapeColor = colour


def lead_colours(map_path=None):
    """{tissue: [r, g, b]} for the two lead tissues, read from tissue_map.yaml.

    Read-only, and via the repo's one parser, so a preview cannot show a colour
    that disagrees with what apply_colors.py would paint or what Ansys is told
    the material is. Returns {} rather than raising if the map is unreadable --
    a preview is worth having in grey.
    """
    try:
        ansys = os.path.join(REPO, "src", "ansys")
        if ansys not in sys.path:
            sys.path.insert(0, ansys)
        from check_tissue_map import parse_tissue_map
        tissues = parse_tissue_map(map_path or os.path.join(ansys, "tissue_map.yaml"))
        return {t["name"]: t["color_rgb"] for t in tissues
                if t["name"] in ("electrode_contact", "lead_insulation")
                and t["color_rgb"]}
    except Exception:
        return {}


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--configs", default=DEFAULT_CONFIGS)
    ap.add_argument("--corridor", default=DEFAULT_CORRIDOR)
    ap.add_argument("--stl-dir", default=DEFAULT_STL)
    ap.add_argument("--out", default=DEFAULT_OUT)
    ap.add_argument("--report", default=DEFAULT_REPORT)
    ap.add_argument("--name", help="a configuration from lead_configs.yaml")
    ap.add_argument("--all", action="store_true", help="every configuration in the file")
    ap.add_argument("--list", action="store_true", help="list them; no FreeCAD needed")
    ap.add_argument("--show", action="store_true",
                    help="resolve to numbers and stop; no FreeCAD needed")
    ap.add_argument("--validate", action="store_true", help="measure the fit")
    ap.add_argument("--export", action="store_true",
                    help="validate, then write STLs if it passes")
    ap.add_argument("--force", action="store_true",
                    help="export even if validation failed (it will not fit)")
    ap.add_argument("--tag", help="output directory name; defaults to the config name")
    ap.add_argument("--max-samples", type=int, default=6000)
    for flag, kind in (("--side", str), ("--contacts", int), ("--contact-length", float),
                       ("--gap", float), ("--diameter", float), ("--x-offset", float),
                       ("--z-center", float), ("--tail", float), ("--level", str)):
        ap.add_argument(flag, type=kind, default=None, help="override the configuration")
    args = ap.parse_args(script_args())

    out = []
    def say(fmt, *a):
        out.append(fmt % a if a else fmt)

    cfg = parse_config_file(args.configs)
    overrides = {k: getattr(args, k.replace("-", "_"))
                 for k in ("side", "contacts", "contact_length", "gap", "diameter",
                           "x_offset", "z_center", "tail", "level")}

    if args.list:
        say("configurations in %s", os.path.relpath(args.configs, REPO))
        say("")
        for name in cfg["configs"]:
            params = resolve(name, cfg)
            say("%-24s %s", name, describe(params))
        say("")
        say("vertebral levels (z of the body centre, mm)")
        for level, z in cfg["levels"].items():
            say("   %-14s %8.2f", level, float(z))
        say("")
        say("RESULT: %d configurations", len(cfg["configs"]))
        write_report(args.report, out)
        return 0

    names = list(cfg["configs"]) if args.all else [args.name]
    if names == [None]:
        say("Nothing to do: give --name NAME, or --all, or --list.")
        write_report(args.report, out)
        return 2

    # Fail here rather than three functions deep with an ImportError traceback:
    # measuring anything needs FreeCAD, and the fix is a different interpreter,
    # not a different configuration.
    if args.validate or args.export:
        try:
            _geometry_modules()
        except ImportError as exc:
            say("Measuring a lead needs FreeCAD, and this interpreter has none (%s).",
                exc)
            say("")
            say("Run it under freecadcmd instead:")
            say('    BUILD_LEAD_CONFIG_ARGS="%s" freecadcmd src/freecad/%s',
                " ".join(script_args()), os.path.basename(__file__))
            say("")
            say("--list and --show need no FreeCAD and work here.")
            say("")
            say("RESULT: aborted, no FreeCAD in this interpreter")
            write_report(args.report, out)
            return 2

    rc = 0
    for name in names:
        try:
            params = resolve(name, cfg, overrides)
        except (KeyError, ValueError) as exc:
            # A bad name or an impossible lead is a user mistake, not a bug:
            # say what is wrong and move on to the next configuration.
            say("=" * 72)
            say("%s", name or "custom")
            say("   cannot use this configuration: %s", exc)
            say("")
            rc = 2
            continue
        say("=" * 72)
        say("%s", name or "custom")
        say("%s", describe(params))
        if args.show:
            for key in sorted(PARAM_TYPES):
                if key in params:
                    say("   %-16s %s", key, params[key])
            say("   %-16s %.2f", "x (absolute)", params["x"])
            say("")
            continue

        result = validate(params, meshes=load_meshes(args.stl_dir),
                          corridor=json.load(open(args.corridor)),
                          max_samples=args.max_samples)
        say("")
        for line in report_lines(params, result)[2:]:
            say("%s", line)
        if not result["ok"]:
            rc = 1

        if args.export:
            if not result["ok"] and not args.force:
                say("")
                say("NOT EXPORTED: this configuration does not fit. Fix it, or pass "
                    "--force if you know what you are doing.")
            else:
                outdir, files = export_stls(params, result["shapes"], args.out, args.tag)
                say("")
                say("wrote %d files to %s", len(files), os.path.relpath(outdir, REPO))
        say("")

    say("RESULT: %s", {0: "OK",
                       1: "one or more configurations do not fit",
                       2: "one or more configurations could not be used"}[rc])
    write_report(args.report, out)
    return rc


def write_report(path, lines):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    if not is_freecad_interpreter():
        sys.stdout.write("\n".join(lines) + "\n")


def run_as_freecadcmd_script():
    """True when a FreeCAD interpreter was handed THIS file to run; see apply_colors.py."""
    if len(sys.argv) < 2 or not is_freecad_interpreter():
        return False
    return os.path.abspath(sys.argv[1]) == os.path.abspath(__file__)


if __name__ == "__main__":
    sys.exit(main())
elif run_as_freecadcmd_script():
    # No sys.exit: SystemExit out of an imported module makes freecadcmd print a
    # traceback over a run that succeeded.
    main()
