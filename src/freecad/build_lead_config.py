#!/usr/bin/env python3
"""Build and VALIDATE any named lead configuration from lead_configs.yaml.

WHAT A CLINICIAN ASKS FOR, AND WHAT THIS TURNS IT INTO
------------------------------------------------------
The clinical description of a lead is "8 contacts, 3 mm long, 1 mm apart,
dorsal, at T10, 1 mm left of midline", or "four contacts on the left third
ganglion". None of that is geometry. This module is the layer that turns it into
geometry and then, more importantly, tells you whether the answer is
anatomically possible in THIS model:

    lead_configs.yaml  ->  resolve_set()   ->  a list of leads, each a flat set
                                               of numbers
                       ->  validate_set()  ->  PASS/WARN/FAIL per lead AND for
                                               the set as a whole
                       ->  export_set()    ->  generated_leads/<name>/*.stl
                       ->  preview_set()   ->  disposable bodies in the open
                                               document

src/freecad/lead_designer.FCMacro is the same thing with spinboxes on it.

A CONFIGURATION IS A SET OF LEADS
---------------------------------
Clinicians do not implant one lead and stop: two dorsal leads straddling the
midline is routine, so is a dorsal lead plus a DRG lead, so is one DRG lead per
side. A configuration here is therefore a SET, capped at `max_leads_per_type`
leads of any one type (default 2). A configuration that describes its lead
inline, with no `leads:` list, is a set of one -- which is why every entry that
existed before sets did still resolves to exactly the geometry it always did.

THREE LEAD TYPES, AND DRG IS A DIFFERENT GEOMETRY PROBLEM
---------------------------------------------------------
    dorsal / ventral    swept along the epidural canal centreline y(z) measured
                        at the requested x by measure_corridor.py. The lead runs
                        rostro-caudally; z is the sweep parameter and x is fixed.

    drg                 swept along the FORAMINAL corridor (y(u), z(u)) measured
                        by measure_foramen.py, where u is the lateral distance
                        from the midline. The lead runs OUT THROUGH THE FORAMEN
                        to one named ganglion: its long axis is x, and both of
                        the other coordinates follow a curve. Nothing about the
                        canal centreline applies to it.

The two are not two settings of one thing, so their parameters do not overlap
and resolve_lead() REJECTS rather than ignores one that means nothing for the
type -- `z_center` on a DRG lead, or `target` on a dorsal one, is an error, on
the same principle that a typo'd `diamter: 2.0` is an error rather than a
silently kept default.

WHY LATERAL OFFSETS RE-MEASURE THE CENTRELINE
---------------------------------------------
epidural_corridor.json holds the centreline y(z) measured at the MIDLINE. The
channel centre moves as you go lateral -- about 1.1 mm of y at 4 mm off the
midline -- so a lead swept along the midline curve at x = midline + 4 sits about
half its clearance off centre before it has done anything else. Any configuration
with a non-zero x_offset therefore gets its own centreline, ray cast and fitted
at ITS x over ITS z span. At x_offset = 0 this reproduces the stored polynomial;
the report prints the agreement so the two can never silently diverge.

The DRG side does the same thing for the same reason: the stored quartic in
foraminal_corridors.json is the reference, and the block is re-measured along
the lead's own u span, with the disagreement reported.

WHY THERE IS NO BOOLEAN SUBTRACTION HERE
----------------------------------------
The intuitive mental model is that a lead has to be CARVED out of the tissue
around it -- Mesh.difference() the epidural fat against the lead so the two do
not occupy the same space. That is not what this model needs, for two reasons,
and the fit check below is what replaces it.

1. RADO's compartments are HOLLOW SHELLS that tile space. Each body already has
   the inner ones carved out, and the epidural fat is a real void that the lead
   sits INSIDE -- the same void a percutaneous lead is threaded into in theatre.
   Ray casting the epidural mesh (see measure_corridor.py) measures that void
   directly: 2.27-2.35 mm dorsally, 1.67-1.74 mm ventrally at the midline. A
   1.3 mm lead fits either side with clearance to spare.

2. Mesh.difference() is a COMPLETE NO-OP in this FreeCAD build -- verified. It
   returns without error and changes nothing. Anything built on it would look
   like it worked and quietly produce an unmodified model.

So the correct check is not "did the subtraction succeed" but "is carving
unnecessary". For an epidural lead that is three measurements:

    fit         the channel thickness measured AT the requested side, lateral
                offset and z span, against the lead diameter
    containment every sampled vertex of the built lead tested for being inside
                the epidural fat and outside the dura
    coverage    whether the lead runs off the end of the measured centreline

A DRG lead asks a different question, because the foramen is not a void: RADO
models the tissue filling each neuroforamen as one SOLID block, and a
percutaneous DRG lead is threaded through foraminal fat rather than dropped into
a gap. So the checks become

    fit         the clearance from the lead axis to the wall of the foraminal
                block, in y and in z, at every station along the lead
    containment every sampled vertex inside that block
    clearance   every sampled vertex OUTSIDE the target ganglion, outside every
                nerve-root core near the corridor, and outside the thecal sac --
                a lead that passes through the ganglion it is meant to stimulate
                is not a placement, it is an injury
    coverage    whether the array runs off the measured corridor

AND ONE CHECK THE SET OWNS RATHER THAN ANY LEAD
-----------------------------------------------
Two leads can each pass every check above and still be un-implantable, because
nothing in a per-lead check looks at the OTHER lead. Two 1.30 mm dorsal leads
1.0 mm apart both fit the channel; they simply cannot both be there. So
validate_set() measures every pair -- the exact shape-to-shape distance, and,
when that is zero, the volume they share and how much of one lead's sampled
surface is inside the other. `dorsal_pair_tight_FAILS` in the catalogue is the
worked example, the way `ventral_z110_2mm_FAILS` is for the fit check.

A configuration that fails is REPORTED, in clinical terms, and not exported. It
is not carved into the anatomy to make it fit.

ONE MODEL DOCUMENT
------------------
NBF_RADO-SCS.FCStd holds anatomy and nothing else. There is no .FCStd per
configuration and no longer a .FCStd per study arm: the durable artefacts are
the YAML entry and generated_leads/<name>/, both small and both text-ish, and
the bodies in a FreeCAD document are a DISPOSABLE PREVIEW:

    preview_set(doc, leads)   adds every lead to the open document inside a
                              group named SCS_LeadPreview, replacing any
                              previous preview
    clear_preview(doc)        deletes that group and everything in it

Explore ten configurations, keep the one you want, export its STLs, and close
the document WITHOUT SAVING.

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
    python3 src/freecad/build_lead_config.py --name dorsal_pair_T11 --show

The measuring and building half does:

    BUILD_LEAD_CONFIG_ARGS="--name drg_L3_4c --validate" \\
        freecadcmd src/freecad/build_lead_config.py
    BUILD_LEAD_CONFIG_ARGS="--all --validate" \\
        freecadcmd src/freecad/build_lead_config.py
    BUILD_LEAD_CONFIG_ARGS="--name dorsal_z110_8c --export" \\
        freecadcmd src/freecad/build_lead_config.py
    BUILD_LEAD_CONFIG_ARGS="--name rado_drg_L3_match --compare-rado" \\
        freecadcmd src/freecad/build_lead_config.py

Any parameter can be overridden on top of a named entry, which is how the macro
drives it for a custom lead. Overrides apply to EVERY lead in the set, so they
are for single-lead configurations and for sweeping one number across a pair:

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
DEFAULT_FORAMEN = os.path.join(HERE, "foraminal_corridors.json")
DEFAULT_STL = os.path.join(REPO, "STL_files")
DEFAULT_OUT = os.path.join(HERE, "generated_leads")
DEFAULT_REPORT = os.path.join(HERE, "build_lead_config_report.txt")

MIDLINE_X = 56.60
EPIDURAL_STL = "T8-10 - neuro_EpiduralSpace-1.STL"
DURA_STL = "T8-10 - neuro_Meninges-1.STL"

# RADO's own 4-contact DRG lead, as it ships in STL_files/. Not part of any
# generated configuration -- it is the independent reference the DRG trajectory
# is checked against by compare_rado(). See rado_drg_L3_match in the catalogue.
RADO_DRG_STLS = ["SCS Lead Electrode %d.stl" % i for i in (1, 2, 3, 4)]
RADO_DRG_TARGET = "L3"

# The name of the group every previewed body lives in. Matched on Name, not
# Label, so retitling it in the tree does not orphan the preview.
PREVIEW_GROUP = "SCS_LeadPreview"

# Clearance below which a fit is reported as tight rather than comfortable. Half
# of it is the amount the lead can wander off the channel centre before it
# touches a wall, and the centreline fit itself is good to about 0.02 mm, so
# 0.20 mm is roughly "ten times the modelling error" rather than a clinical rule.
TIGHT_CLEARANCE_MM = 0.20

# The same idea for two leads beside each other: below this much fat between
# them they are reported as touching-tight rather than comfortably apart. Wider
# than TIGHT_CLEARANCE_MM because two leads can each drift by their own
# centreline error, and because a mesher given a 0.2 mm gap between two curved
# conductors produces slivers.
TIGHT_LEAD_GAP_MM = 0.40

# Every parameter a lead may carry, and how to read it. Anything not in here is
# rejected by resolve_lead() rather than silently ignored -- a typo'd
# "diamter: 2.0" that quietly kept the 1.3 mm default would be the worst kind of
# bug here.
PARAM_TYPES = {
    "type": str, "contacts": int, "contact_length": float, "gap": float,
    "diameter": float, "tail": float, "note": str, "label": str,
    # dorsal / ventral only
    "side": str, "x_offset": float, "z_center": float, "level": str,
    # DRG only
    "target": str, "lateral_offset": float, "ganglion_clearance": float,
    "y_offset": float, "z_offset": float,
}

# Which parameters each type accepts. The split is the point: a DRG lead has no
# z_center and a dorsal lead has no target, and saying so is what stops a
# configuration from carrying a number that quietly does nothing.
COMMON_PARAMS = {"type", "contacts", "contact_length", "gap", "diameter",
                 "tail", "note", "label"}
EPIDURAL_PARAMS = {"side", "x_offset", "z_center", "level"}
DRG_PARAMS = {"target", "lateral_offset", "ganglion_clearance",
              "y_offset", "z_offset"}

TYPES = ("dorsal", "ventral", "drg")
EPIDURAL_TYPES = ("dorsal", "ventral")

DEFAULT_MAX_LEADS_PER_TYPE = 2


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
    """Read lead_configs.yaml.

    Returns the top-level mapping: scalars like `max_leads_per_type` as scalars,
    and `defaults`, `epidural_defaults`, `drg_defaults`, `levels` and `configs`
    as mappings. A config entry is a mapping of scalars, optionally with a
    `leads:` key holding a list of mappings.

    A hand-rolled parser, for the reason requirements.txt gives: PyYAML is
    deliberately not a dependency of this directory, so the scripts run on a bare
    interpreter with nothing installed. It accepts exactly the subset
    lead_configs.yaml uses and nothing else -- no flow style, no anchors, no
    multi-line scalars, and lists only in the one place they are needed.
    Anything it cannot parse RAISES rather than being skipped, because a
    silently dropped parameter is a wrong lead.

    The grammar, by indent:

        0   key: value          a top-level scalar
        0   key:                a section
        2   key: value          a scalar in the section
        2   key:                an entry in the section
        4   key: value          a scalar in the entry
        4   leads:              the one list, and only under `configs`
        6   - key: value        the first line of one lead in that list
        8   key: value          a further line of that lead
    """
    top, section, entry, leads, lead = {}, None, None, None, None
    for n, raw in enumerate(open(path), 1):
        line = raw.rstrip("\n")
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        indent = len(line) - len(line.lstrip())

        if indent == 6 and line.lstrip().startswith("- "):
            if leads is None:
                raise ValueError("%s:%d: list item outside a `leads:` list" % (path, n))
            lead = {}
            leads.append(lead)
            body = line.lstrip()[2:]
        elif indent == 8:
            if lead is None:
                raise ValueError("%s:%d: continuation outside a list item" % (path, n))
            body = line.lstrip()
        else:
            body = line.lstrip()

        m = re.match(r"^([A-Za-z_][A-Za-z0-9_.\-]*)\s*:\s*(.*)$", body)
        if not m:
            raise ValueError("%s:%d: cannot parse %r" % (path, n, line))
        key, rest = m.group(1), strip_comment(m.group(2))

        if indent in (6, 8):
            if rest == "":
                raise ValueError("%s:%d: %r inside a lead must have a value"
                                 % (path, n, key))
            lead[key] = scalar(rest)
            continue

        lead = None
        if indent == 0:
            leads = None
            if rest:
                top[key] = scalar(rest)
                section, entry = None, None
            else:
                section, entry = {}, None
                top[key] = section
        elif indent == 2:
            leads = None
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
            if key == "leads":
                if rest:
                    raise ValueError("%s:%d: `leads:` must be a list" % (path, n))
                leads = []
                entry["leads"] = leads
            elif rest == "":
                raise ValueError("%s:%d: %r has no value" % (path, n, key))
            else:
                leads = None
                entry[key] = scalar(rest)
        else:
            raise ValueError("%s:%d: unexpected indent %d" % (path, n, indent))

    for required in ("defaults", "levels", "configs"):
        if required not in top:
            raise ValueError("%s: no %r section" % (path, required))
    top.setdefault("epidural_defaults", {})
    top.setdefault("drg_defaults", {})
    top.setdefault("max_leads_per_type", DEFAULT_MAX_LEADS_PER_TYPE)
    return top


def strip_comment(text):
    """Drop a trailing # comment, unless it is inside quotes."""
    out, quote = [], None
    for ch in text:
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


# --------------------------------------------------------------------------
# resolving a configuration into leads
# --------------------------------------------------------------------------
def config_names(cfg):
    return list(cfg["configs"])


def lead_entries(cfg, name):
    """The per-lead mappings of one configuration, and its shared keys.

    Returns (shared, [lead, ...]). A configuration with no `leads:` list is a set
    of one whose single lead is the entry itself -- which is what makes every
    pre-existing single-lead entry resolve to exactly the geometry it always did.
    """
    if name is None:
        return {}, [{}]
    if name not in cfg["configs"]:
        raise ValueError("no configuration named %r (have: %s)"
                         % (name, ", ".join(sorted(cfg["configs"]))))
    entry = dict(cfg["configs"][name])
    leads = entry.pop("leads", None)
    if leads is None:
        return {}, [entry]
    if not leads:
        raise ValueError("configuration %r has an empty `leads:` list" % name)
    return entry, [dict(l) for l in leads]


def resolve_set(name, cfg, overrides=None):
    """Flatten one named configuration into a list of leads.

    Each lead is resolved independently through resolve_lead(); the set-level
    rule enforced here is the only thing a lead cannot know on its own, which is
    how many leads of its type the configuration is allowed to hold.
    """
    shared, entries = lead_entries(cfg, name)
    leads = []
    for index, lead_entry in enumerate(entries, start=1):
        leads.append(resolve_lead(cfg, name, shared, lead_entry, overrides, index,
                                  len(entries)))

    cap = int(cfg.get("max_leads_per_type", DEFAULT_MAX_LEADS_PER_TYPE))
    counts = {}
    for lead in leads:
        counts[lead["type"]] = counts.get(lead["type"], 0) + 1
    over = sorted(t for t, c in counts.items() if c > cap)
    if over:
        raise ValueError(
            "configuration %r asks for %s, but max_leads_per_type is %d"
            % (name or "custom",
               " and ".join("%d %s leads" % (counts[t], t) for t in over), cap))
    return leads


def resolve_lead(cfg, name, shared, entry, overrides, index=1, of=1):
    """Flatten ONE lead into the numbers the geometry needs.

    Layering, later wins:

        defaults  ->  the type's own default block  ->  the configuration's
        shared keys  ->  this lead's own keys  ->  overrides (CLI flags, or the
        macro's spinboxes)

    The type is settled first, across every layer, because it decides which
    default block sits at position two and which parameters are even legal.

    `level: T10` is turned into the z of that vertebral body here, so nothing
    downstream has to know what a level is; an explicit z_center in a LATER layer
    beats an earlier level, and vice versa, which is what lets the macro offer
    both controls without them fighting.

    Returns a plain dict, plus "name", "index", "label", and for an epidural lead
    "x" (absolute), for a DRG lead the ganglion "target".
    """
    ovr = {k: v for k, v in (overrides or {}).items() if v is not None}

    kind = lead_type(cfg["defaults"], shared, entry, ovr)
    type_defaults = cfg["drg_defaults"] if kind == "drg" else cfg["epidural_defaults"]
    layers = [cfg["defaults"], type_defaults, shared, entry, ovr]

    params = {}
    for layer in layers:
        params.update(layer)
    params.pop("side", None)
    params["type"] = kind

    allowed = COMMON_PARAMS | (DRG_PARAMS if kind == "drg" else EPIDURAL_PARAMS)
    unknown = [k for k in params if k not in PARAM_TYPES]
    if unknown:
        raise ValueError("unknown parameter(s) %s in configuration %r"
                         % (", ".join(sorted(unknown)), name or "custom"))
    wrong = [k for k in params if k not in allowed]
    if wrong:
        other = "a DRG lead" if kind != "drg" else "a dorsal or ventral lead"
        raise ValueError(
            "%s is a %s lead, so %s mean%s nothing to it -- %s belong%s to %s"
            % (name or "custom", kind, ", ".join(sorted(wrong)),
               "" if len(wrong) > 1 else "s", "they" if len(wrong) > 1 else "it",
               "" if len(wrong) > 1 else "s", other))

    for k, v in list(params.items()):
        if isinstance(v, str) and PARAM_TYPES[k] is not str:
            params[k] = PARAM_TYPES[k](v)

    if kind in EPIDURAL_TYPES:
        level = params.get("level")
        if level:
            if level not in cfg["levels"]:
                raise ValueError("no vertebral level named %r (have: %s)"
                                 % (level, ", ".join(cfg["levels"])))
            # A level named in the same or a later layer than an explicit
            # z_center wins; that is what "later wins" means when both are
            # present, and the only way to tell is which layer each came from.
            if _layer_of("z_center", layers) <= _layer_of("level", layers):
                params["z_center"] = float(cfg["levels"][level])
        params["x"] = MIDLINE_X + float(params["x_offset"])

    params["name"] = name or "custom"
    params["index"] = index
    params["of"] = of
    params.setdefault("label", default_label(params, index, of))
    check_sane(params)
    return params


def lead_type(*layers):
    """The lead's type, from the last layer that names one. Accepts `side:`.

    `side: dorsal | ventral` was what the catalogue said before DRG leads
    existed, and both the CLI and older entries may still use it, so it is
    accepted as a spelling of `type` rather than becoming a confusing "unknown
    parameter" error. It cannot spell `drg`: a DRG lead is not a side.
    """
    kind = "dorsal"
    for layer in layers:
        if layer.get("type") is not None:
            kind = str(layer["type"])
        if layer.get("side") is not None:
            side = str(layer["side"])
            if side not in EPIDURAL_TYPES:
                raise ValueError(
                    "side must be dorsal or ventral, not %r -- a DRG lead is not "
                    "a side of the canal, write `type: drg`" % side)
            kind = side
    if kind not in TYPES:
        raise ValueError("type must be one of %s, not %r"
                         % (", ".join(TYPES), kind))
    return kind


def default_label(params, index, of):
    """A short human name for one lead of a set, used in the tree and on disk."""
    if of == 1:
        return "%s lead" % params["type"]
    return "lead %d (%s)" % (index, params["type"])


def _layer_of(key, layers):
    """Index of the last layer that set `key`; -1 if none did. See resolve_lead()."""
    last = -1
    for i, layer in enumerate(layers):
        if layer.get(key) is not None:
            last = i
    return last


def check_sane(params):
    """Reject nonsense before any geometry is attempted, in the user's terms."""
    problems = []
    if params["contacts"] < 1:
        problems.append("a lead needs at least one contact")
    for key, what in (("contact_length", "contact length"), ("diameter", "diameter")):
        if params[key] <= 0:
            problems.append("%s must be greater than zero" % what)
    if params["gap"] < 0:
        problems.append("the gap between contacts cannot be negative")
    if params["tail"] < 0:
        problems.append("the tail cannot be negative")
    if params["type"] == "drg":
        target = str(params.get("target", ""))
        if not re.match(r"^[LR][1-9]\d*$", target):
            problems.append("target must name a ganglion, like L3 or R2, not %r"
                            % target)
    if problems:
        raise ValueError("; ".join(problems))


def array_length(params):
    """Length of the contact array, mm (first contact's start to last one's end)."""
    return (params["contacts"] * params["contact_length"]
            + (params["contacts"] - 1) * params["gap"])


def total_length(params):
    return array_length(params) + 2 * params["tail"]


def describe(params):
    """One lead as a clinician would say it out loud."""
    kind = params["type"]
    common = ("%d contacts, %.2f mm long, %.2f mm apart, %.2f mm diameter"
              % (params["contacts"], params["contact_length"], params["gap"],
                 params["diameter"]))
    if kind == "drg":
        target = params["target"]
        side = "left" if target[0] == "L" else "right"
        shift = float(params.get("lateral_offset", 0.0))
        if abs(shift) < 1e-9:
            where = "centred on the ganglion"
        else:
            where = "%.1f mm %s the ganglion along the root" % (
                abs(shift), "distal to" if shift > 0 else "proximal to")
        return ("%s, threaded through the %s neuroforamen to ganglion %s, %s"
                % (common, side, target, where))
    dx = float(params["x_offset"])
    if abs(dx) < 1e-9:
        where = "on the midline"
    else:
        where = "%.1f mm %s of midline" % (abs(dx), "left" if dx > 0 else "right")
    at = params.get("level") or "z = %.1f mm" % params["z_center"]
    return ("%s, in the %s epidural space at %s, %s" % (common, kind, at, where))


def describe_set(leads):
    """The whole configuration in one line."""
    if len(leads) == 1:
        return describe(leads[0])
    counts = {}
    for lead in leads:
        counts[lead["type"]] = counts.get(lead["type"], 0) + 1
    return "%d leads: %s" % (
        len(leads), ", ".join("%d %s" % (counts[t], t) for t in TYPES if t in counts))


# --------------------------------------------------------------------------
# measurement and geometry -- everything below needs FreeCAD
# --------------------------------------------------------------------------
def _geometry_modules():
    """Import the FreeCAD-only modules. Deferred so --list/--show run anywhere.

    make_scs_lead, measure_corridor and measure_foramen all import FreeCAD at
    module scope, and all call main() at module scope UNLESS handed to freecadcmd
    as the script -- which is why they each carry a run_as_freecadcmd_script()
    guard. Without it, importing measure_corridor here would silently re-measure
    the corridor and rewrite epidural_corridor.json.
    """
    if HERE not in sys.path:
        sys.path.insert(0, HERE)
    import FreeCAD
    import Mesh
    import Part
    import make_scs_lead
    import measure_corridor
    import measure_foramen
    return FreeCAD, Mesh, Part, make_scs_lead, measure_corridor, measure_foramen


def load_meshes(stl_dir=DEFAULT_STL):
    """The two meshes every epidural check needs: the epidural fat and the dura.

    Returned as a plain dict which the DRG path then adds foraminal blocks and
    neural bodies to, keyed by filename, so one caller-held cache serves a set
    containing leads of all three types.
    """
    _fc, Mesh, _p, _msl, _mc, _mf = _geometry_modules()
    return {"epidural": Mesh.Mesh(os.path.join(stl_dir, EPIDURAL_STL)),
            "dura": Mesh.Mesh(os.path.join(stl_dir, DURA_STL))}


def mesh_named(meshes, filename, stl_dir=DEFAULT_STL):
    """One STL, loaded once per `meshes` dict and kept."""
    _fc, Mesh, _p, _msl, _mc, _mf = _geometry_modules()
    if filename not in meshes:
        meshes[filename] = Mesh.Mesh(os.path.join(stl_dir, filename))
    return meshes[filename]


def load_foramen_corridors(path=DEFAULT_FORAMEN):
    """foraminal_corridors.json, with a useful error when it has not been made."""
    if not os.path.exists(path):
        raise ValueError(
            "no foraminal corridor measurements at %s -- a DRG lead needs them. "
            "Run: freecadcmd src/freecad/measure_foramen.py"
            % os.path.relpath(path, REPO))
    return json.load(open(path))


def corridor_for(params, foramen):
    """The measured corridor for this lead's target ganglion."""
    target = params["target"]
    if target not in foramen["ganglia"]:
        raise ValueError("no ganglion named %r in this model (have: %s)"
                         % (target, ", ".join(sorted(foramen["ganglia"]))))
    return foramen["ganglia"][target]


# -- epidural (dorsal / ventral) -------------------------------------------
def channel_at(epi, x, z, side):
    """The epidural channel on one side at (x, z): (y_lo, y_hi, thickness).

    None when the ray does not give exactly two material intervals. Two is the
    signature of a clean slice -- ventral fat, then the thecal sac, then dorsal
    fat -- and anything else means the ray went through a foramen, a root sleeve
    or off the end of the sac, where "the thickness of the channel" is not a
    well-defined quantity. Reported rather than averaged over.
    """
    _fc, _Mesh, _p, _msl, mc, _mf = _geometry_modules()
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
    _fc, _Mesh, _p, _msl, mc, _mf = _geometry_modules()
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
    _fc, _Mesh, _p, _msl, mc, _mf = _geometry_modules()
    side = params["type"]
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


# -- foraminal (DRG) --------------------------------------------------------
def drg_point_of(params, corridor, standoff=0.0):
    """The curve a DRG lead is swept along: u -> (x, y, z).

    The stored quartic from foraminal_corridors.json, shifted by three things:
    the derived `standoff` (see ganglion_standoff), and the lead's own y_offset
    and z_offset. The last two are how a lead is nudged off the corridor by hand
    -- dorsal/ventral and rostral/caudal -- and they are the DRG analogue of
    x_offset on an epidural lead. Unlike x_offset they do NOT change which curve
    is measured: the corridor is a property of the foramen, and sliding a lead
    within it does not move its walls.
    """
    _fc, _Mesh, _p, _msl, _mc, mf = _geometry_modules()
    cy = corridor["y_poly_coef_high_to_low"]
    cz = corridor["z_poly_coef_high_to_low"]
    sign = corridor["x_sign"]
    dy = float(params.get("y_offset", 0.0))
    dz = float(params.get("z_offset", 0.0)) + standoff
    return lambda u: (MIDLINE_X + sign * u,
                      mf.polyval(cy, u) + dy,
                      mf.polyval(cz, u) + dz)


def ganglion_standoff(params, corridor, drg, u_lo, u_hi, step=0.5):
    """How far ROSTRALLY the lead must ride to clear its ganglion. mm, >= 0.

    WHY THIS IS DERIVED AND NOT A NUMBER IN THE YAML
    ------------------------------------------------
    The centre of the foraminal corridor runs almost exactly along the TOP of
    the ganglion -- measured on L3 it passes within about 0.1 mm of the ganglion
    surface, and RADO's own hand-placed lead sits in the same fraction of a
    millimetre. That is not an accident and it is not a problem: a DRG lead is
    supposed to lie against its ganglion, because the whole point is a low
    activation threshold. But it means a lead dropped exactly on the corridor
    centre is a coin toss between grazing the ganglion and cutting into it.

    So the clinical parameter is not "shift the lead 1.4 mm up", it is "hold the
    lead this far off the ganglion", and the shift is solved for. It has to be,
    because the answer differs by foramen: holding 0.25 mm of clearance needs
    +0.77 mm on L3 and +1.62 mm on R3, and either number written into the
    catalogue would be wrong for the other side.

    ROSTRALLY (+z, towards the pedicle above) rather than dorsally, for two
    reasons. It is where the room is -- the foraminal block is 13-14 mm tall in
    z against 8-10 mm in y, and the ganglion sits low in it, leaving 2.5 mm of
    headroom above the lead and much less to either side. And it is where a
    percutaneous DRG lead goes in theatre: into the superior aspect of the
    foramen, under the pedicle, over the ganglion.

    Measured by casting +Z through the ganglion at three y across the lead's own
    width, so a ganglion whose top is higher at the edge of the lead than under
    its axis is still cleared.
    """
    _fc, _Mesh, _p, _msl, _mc, mf = _geometry_modules()
    want = float(params.get("ganglion_clearance", 0.0))
    if want < 0:
        raise ValueError("ganglion_clearance cannot be negative: %r" % want)
    cy = corridor["y_poly_coef_high_to_low"]
    cz = corridor["z_poly_coef_high_to_low"]
    sign = corridor["x_sign"]
    radius = 0.5 * float(params["diameter"])
    dy = float(params.get("y_offset", 0.0))

    need = 0.0
    u = u_lo
    while u <= u_hi + 1e-9:
        x = MIDLINE_X + sign * u
        y0 = mf.polyval(cy, u) + dy
        z0 = mf.polyval(cz, u) + float(params.get("z_offset", 0.0))
        for y in (y0 - radius, y0, y0 + radius):
            segs = mf.ray_intervals(drg, (x, y, 0.0), (0.0, 0.0, 1.0))
            if segs:
                top = max(hi for _lo, hi in segs)
                need = max(need, top + want + radius - z0)
        u += step
    return max(0.0, need)


def drg_u_center(params, corridor):
    """Where the centre of the contact array sits along the corridor.

    Zero lateral_offset centres the array on the ganglion, which is the
    placement a DRG lead is aiming for: the contacts straddle the ganglion so
    the field covers it. Positive is further out along the root.
    """
    return float(corridor["ganglion_u_centre"]) + float(params.get("lateral_offset", 0.0))


def remeasure_foramen(block, corridor, params, u_lo, u_hi, step=0.5):
    """Walk the foraminal block along the lead's own span. Returns (stations, info).

    stations is [(u, y_centre, z_centre, y_extent, z_extent)] and info carries
    the worst disagreement with the stored quartic, which is the live check that
    the stored fit and the block have not drifted apart -- the same role the
    midline comparison plays for an epidural lead.
    """
    _fc, _Mesh, _p, _msl, _mc, mf = _geometry_modules()
    cy = corridor["y_poly_coef_high_to_low"]
    cz = corridor["z_poly_coef_high_to_low"]
    _gx, gy, gz = corridor["ganglion_centroid"]
    rows = mf.scan_corridor(block, corridor["x_sign"], gy, gz, u_lo, u_hi, step)
    if not rows:
        raise ValueError(
            "the %s foraminal block could not be measured over u %.1f..%.1f -- "
            "the lead is off the end of the foramen"
            % (params["target"], u_lo, u_hi))
    # Compared only over the range the quartic was FITTED to. Outside it the
    # polynomial is extrapolating into the block's taper, where the measured
    # centre swings by a millimetre or more over a millimetre of u -- a real
    # property of the taper, not a drift between the fit and the block, and
    # folding it into this number would make the check cry wolf on every lead
    # whose tip runs a little past the ganglion.
    lo, hi = corridor["usable_u"]
    fitted = [r for r in rows if lo - 1e-9 <= r[0] <= hi + 1e-9]
    diff = max((max(abs(r[1] - mf.polyval(cy, r[0])),
                    abs(r[2] - mf.polyval(cz, r[0]))) for r in fitted), default=0.0)
    return rows, {"max_disagreement_mm": diff, "n_stations": len(rows),
                  "n_extrapolated": len(rows) - len(fitted)}


def wall_clearance(block, corridor, point_of, u_lo, u_hi, step=0.5):
    """Distance from the lead axis to the foraminal block's wall, station by station.

    Returns [(u, clearance_y, clearance_z)] where each clearance is the SMALLER
    of the two directions -- how far the lead axis can move before it leaves the
    block. This is what "does it fit the foramen" means here, and it is not the
    same as the block's cross-section: a lead pushed to one side of a 7 mm
    channel has 0.5 mm of room even though the channel is roomy.
    """
    _fc, _Mesh, _p, _msl, _mc, mf = _geometry_modules()
    rows = []
    u = u_lo
    while u <= u_hi + 1e-9:
        _x, y, z = point_of(u)
        x = MIDLINE_X + corridor["x_sign"] * u
        sy = mf.containing(mf.ray_intervals(block, (x, 0.0, z), (0.0, 1.0, 0.0)), y)
        sz = mf.containing(mf.ray_intervals(block, (x, y, 0.0), (0.0, 0.0, 1.0)), z)
        if sy is None or sz is None:
            rows.append((u, -1.0, -1.0))
        else:
            rows.append((u, min(y - sy[0], sy[1] - y), min(z - sz[0], sz[1] - z)))
        u += step
    return rows


# --------------------------------------------------------------------------
# building the geometry
# --------------------------------------------------------------------------
def build_shapes(params, path_info):
    """The lead itself: [(kind, index, shape)], plus the segment plan.

    One line of geometry, and it is not written here -- make_scs_lead owns the
    sweep, so all three lead types differ only in the curve handed to it.
    path_info is (point_of, t_center) from whichever measurement applies.
    """
    _fc, _Mesh, _p, msl, _mc, _mf = _geometry_modules()
    point_of, t_center = path_info
    segments, t0, total = msl.segment_plan(
        params["contacts"], params["contact_length"], params["gap"],
        params["tail"], t_center)
    shapes = msl.sweep_path(segments, point_of, 0.5 * params["diameter"])
    return shapes, segments, t0, total


def sample_points(shapes, tolerance=0.02, max_samples=6000):
    """Vertices of the built lead, for the containment and collision checks.

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


def inside_any(mesh, points):
    """The subset of `points` inside `mesh`, bounding-box filtered first.

    The filter is not an optimisation detail: a DRG lead is tested against a
    couple of dozen nerve-root bodies, and almost every point misses almost
    every body's bounding box, so without it the check costs a ray cast per
    point per body.
    """
    bb = mesh.BoundBox
    out = []
    for p in points:
        if not (bb.XMin <= p[0] <= bb.XMax and bb.YMin <= p[1] <= bb.YMax
                and bb.ZMin <= p[2] <= bb.ZMax):
            continue
        if inside(mesh, p):
            out.append(p)
    return out


def min_distance_to_mesh(points, mesh, max_mesh_points=1500):
    """Approximate nearest approach from a point set to a mesh surface, mm.

    Vertex to vertex, with the mesh subsampled by a fixed stride. It is an
    approximation and reads slightly LONG (the true nearest point can be in the
    middle of a facet), so it is reported as a gauge -- "about 1.8 mm clear of
    the ganglion" -- and never as the verdict. The verdict is the parity test in
    inside_any(), which is exact.
    """
    verts = mesh.Points
    stride = max(1, len(verts) // max_mesh_points)
    targets = [(v.x, v.y, v.z) for v in verts[::stride]]
    best = float("inf")
    for p in points:
        for q in targets:
            d = (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2
            if d < best:
                best = d
    return best ** 0.5


# --------------------------------------------------------------------------
# validation -- one lead
# --------------------------------------------------------------------------
def validate_lead(params, meshes=None, corridor=None, foramen=None,
                  max_samples=6000, stl_dir=DEFAULT_STL):
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
        meshes = load_meshes(stl_dir)
    if params["type"] == "drg":
        return _validate_drg(params, meshes, foramen, max_samples, stl_dir)
    if corridor is None:
        corridor = json.load(open(DEFAULT_CORRIDOR))
    return _validate_epidural(params, meshes, corridor, max_samples)


def _validate_epidural(params, meshes, corridor, max_samples):
    epi, dura = meshes["epidural"], meshes["dura"]
    side = params["type"]
    checks, numbers = [], {}

    length = total_length(params)
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
    shapes, segments, _z0, _total = build_shapes(
        params, (lambda z: (params["x"], y_of_z(z), z), params["z_center"]))

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
        checks.append(("FAIL", "Overlaps the dura",
                       "%.1f%% of %d sampled points on the lead surface are inside "
                       "the dura, between z %.1f and %.1f mm. The lead is not in the "
                       "epidural space there -- it is in the thecal sac."
                       % (pct_dura, len(pts), min(p[2] for p in in_dura),
                          max(p[2] for p in in_dura))))
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
    return {"ok": ok, "checks": checks, "numbers": numbers, "points": pts,
            "shapes": shapes, "segments": segments, "params": params}


def _validate_drg(params, meshes, foramen, max_samples, stl_dir):
    """The DRG lead's four checks. See the module docstring for why they differ."""
    checks, numbers = [], {}
    if foramen is None:
        foramen = load_foramen_corridors()
    try:
        corridor = corridor_for(params, foramen)
    except ValueError as exc:
        checks.append(("FAIL", "No such ganglion", str(exc)))
        return {"ok": False, "checks": checks, "numbers": numbers,
                "shapes": None, "params": params}

    target = params["target"]
    u_center = drg_u_center(params, corridor)
    length = total_length(params)
    u0 = u_center - 0.5 * length
    u1 = u0 + length
    c0 = u_center - 0.5 * array_length(params)
    c1 = c0 + array_length(params)
    usable = corridor["usable_u"]
    ganglion_u = corridor["ganglion_u"]
    numbers.update(lead_u=[u0, u1], contacts_u=[c0, c1], length_mm=length,
                   array_mm=array_length(params), usable_u=usable,
                   ganglion_u=ganglion_u, target=target,
                   ganglion_z=corridor["ganglion_centroid"][2],
                   side=corridor["side"])

    # -- 1. does the measured corridor reach this far? --------------------
    #
    # This one only ever WARNS, and that is deliberate. usable_u is the range
    # over which the FITTED centreline is trustworthy, not the extent of the
    # tissue: outside it the quartic is being extrapolated, which is worth
    # saying, but whether the lead is actually somewhere it can be is decided by
    # the containment and wall checks below, which measure the anatomy directly.
    # RADO's own lead is the case that settles it -- its outermost contact sits
    # 26.0 mm from the midline, past L3's measurable channel at 24.5 mm, and it
    # is nevertheless 100% inside the foraminal tissue.
    if c0 < usable[0] or c1 > usable[1]:
        checks.append(("WARN", "Contacts run past the measured corridor",
                       "The contacts sit %.1f-%.1f mm from the midline, and the %s "
                       "foramen is only a cleanly measurable channel over "
                       "%.1f-%.1f mm, so %.1f mm of the active array follows an "
                       "extrapolated curve. The containment check below is what "
                       "says whether it is still inside the foramen."
                       % (c0, c1, target, usable[0], usable[1],
                          max(usable[0] - c0, c1 - usable[1]))))
    elif u0 < usable[0] or u1 > usable[1]:
        checks.append(("WARN", "Tail runs off the measured corridor",
                       "Every contact is inside the measured foramen (%.1f-%.1f mm "
                       "from the midline), but the passive tail reaches %.1f-%.1f mm, "
                       "past the measured range %.1f-%.1f mm, so the last %.1f mm of "
                       "insulator follows an extrapolated curve."
                       % (c0, c1, u0, u1, usable[0], usable[1],
                          max(usable[0] - u0, u1 - usable[1]))))
    else:
        checks.append(("PASS", "Inside the measured foramen",
                       "The lead spans %.1f-%.1f mm from the midline, within the "
                       "%s foramen's measured channel %.1f-%.1f mm."
                       % (u0, u1, target, usable[0], usable[1])))

    # -- 2. build it along the foraminal corridor -------------------------
    block = mesh_named(meshes, corridor["block_stl"], stl_dir)
    drg = mesh_named(meshes, corridor["drg_core_stl"], stl_dir)
    standoff = ganglion_standoff(params, corridor, drg, u0, u1)
    numbers["ganglion_standoff_mm"] = standoff
    numbers["ganglion_clearance_asked_mm"] = float(params.get("ganglion_clearance", 0.0))
    point_of = drg_point_of(params, corridor, standoff)
    shapes, segments, _u0, _total = build_shapes(params, (point_of, u_center))

    # Measured over the LEAD's own span, not the corridor's: clamping to
    # usable_u would hide exactly the overrun the check above only warns about.
    span_lo, span_hi = u0, u1
    try:
        _rows, finfo = remeasure_foramen(block, corridor, params, span_lo, span_hi)
    except ValueError as exc:
        checks.append(("FAIL", "The foramen could not be followed here", str(exc)))
        return {"ok": False, "checks": checks, "numbers": numbers,
                "shapes": shapes, "segments": segments, "params": params}
    numbers["centreline_source"] = ("stored quartic for %s (foraminal_corridors.json, "
                                    "y rms %.3f mm, z rms %.3f mm)"
                                    % (target, corridor["y_fit_rms_mm"],
                                       corridor["z_fit_rms_mm"]))
    numbers["centreline_vs_stored_max_mm"] = finfo["max_disagreement_mm"]
    numbers["extrapolated_stations"] = finfo["n_extrapolated"]

    # -- 3. how much room is there around the lead, along its length? -----
    walls = wall_clearance(block, corridor, point_of, span_lo, span_hi)
    good = [w for w in walls if w[1] >= 0]
    numbers["stations_off_the_block"] = len(walls) - len(good)
    if not good:
        checks.append(("FAIL", "Outside the foramen",
                       "No station along this lead lies inside the %s foraminal "
                       "tissue at all. The offsets have pushed it out of the "
                       "foramen entirely." % target))
        return {"ok": False, "checks": checks, "numbers": numbers, "points": None,
                "shapes": shapes, "segments": segments, "params": params}
    if len(good) < len(walls):
        off = [w[0] for w in walls if w[1] < 0]
        checks.append(("FAIL", "Part of the lead is outside the foramen",
                       "%d of %d stations along the lead axis are outside the %s "
                       "foraminal tissue altogether (%.1f-%.1f mm from the "
                       "midline). That part of the lead is in bone or beyond the "
                       "end of the modelled root."
                       % (len(walls) - len(good), len(walls), target,
                          min(off), max(off))))
    worst = min(good, key=lambda w: min(w[1], w[2]))
    room = min(worst[1], worst[2])
    radius = 0.5 * params["diameter"]
    clearance = room - radius
    numbers.update(wall_room_mm=room, wall_room_at_u=worst[0],
                   clearance_mm=clearance,
                   channel_min_mm=2.0 * room, channel_min_at_u=worst[0])
    if clearance < 0:
        checks.append(("FAIL", "Presses through the wall of the foramen",
                       "At %.1f mm from the midline the lead axis is only %.2f mm "
                       "from the edge of the %s foraminal tissue, and the lead is "
                       "%.2f mm across, so %.2f mm of it is outside. Reduce the "
                       "offsets, or move the array along the root."
                       % (worst[0], room, target, params["diameter"], -clearance)))
    elif clearance < TIGHT_CLEARANCE_MM:
        checks.append(("WARN", "A tight fit in the foramen",
                       "At %.1f mm from the midline the lead has only %.2f mm "
                       "between its surface and the edge of the foraminal tissue."
                       % (worst[0], clearance)))
    else:
        checks.append(("PASS", "Sits clear inside the foramen",
                       "The tightest point is %.1f mm from the midline, where the "
                       "lead axis is %.2f mm from the wall against a %.2f mm lead: "
                       "%.2f mm of clearance."
                       % (worst[0], room, params["diameter"], clearance)))

    # -- 4. containment, and the structures it must not be in -------------
    pts = sample_points(shapes, max_samples=max_samples)
    numbers["sampled_vertices"] = len(pts)
    in_block = inside_any(block, pts)
    pct_block = 100.0 * len(in_block) / len(pts)
    numbers["pct_in_foraminal_tissue"] = pct_block

    dura = meshes["dura"]
    in_dura = inside_any(dura, pts)
    in_drg = inside_any(drg, pts)
    numbers["pct_in_dura"] = 100.0 * len(in_dura) / len(pts)
    numbers["pct_in_ganglion"] = 100.0 * len(in_drg) / len(pts)

    lead_bb = (min(p[0] for p in pts), max(p[0] for p in pts),
               min(p[1] for p in pts), max(p[1] for p in pts),
               min(p[2] for p in pts), max(p[2] for p in pts))
    pierced = []
    for fn in corridor.get("neural_stls", []):
        if fn == corridor["drg_core_stl"]:
            continue
        nerve = mesh_named(meshes, fn, stl_dir)
        bb = nerve.BoundBox
        if (bb.XMin > lead_bb[1] or bb.XMax < lead_bb[0]
                or bb.YMin > lead_bb[3] or bb.YMax < lead_bb[2]
                or bb.ZMin > lead_bb[5] or bb.ZMax < lead_bb[4]):
            continue
        hit = inside_any(nerve, pts)
        if hit:
            pierced.append((fn, len(hit)))
    numbers["nerve_bodies_pierced"] = len(pierced)

    gap_drg = min_distance_to_mesh(pts[::max(1, len(pts) // 400)], drg)
    numbers["gap_to_ganglion_mm"] = gap_drg

    if in_drg:
        checks.append(("FAIL", "Passes through the ganglion",
                       "%.1f%% of %d sampled points on the lead surface are inside "
                       "ganglion %s itself. A DRG lead sits BESIDE its ganglion, in "
                       "the foraminal fat; one that goes through it is not a "
                       "placement. Move it off the corridor centre with z_offset "
                       "(rostral, towards the pedicle) or y_offset."
                       % (numbers["pct_in_ganglion"], len(pts), target)))
    elif in_dura:
        checks.append(("FAIL", "Enters the thecal sac",
                       "%.1f%% of %d sampled points are inside the dura. The lead is "
                       "too medial -- it has been pushed back into the canal rather "
                       "than out through the foramen."
                       % (numbers["pct_in_dura"], len(pts))))
    elif pct_block < 100.0:
        stray = [p for p in pts if p not in in_block]
        checks.append(("FAIL", "Leaves the foraminal tissue",
                       "%.1f%% of %d sampled points on the lead surface are inside "
                       "the %s foraminal tissue; the other %.1f%% are outside it, "
                       "between %.1f and %.1f mm from the midline. That part of the "
                       "lead would be in bone."
                       % (pct_block, len(pts), target, 100.0 - pct_block,
                          min(abs(p[0] - MIDLINE_X) for p in stray),
                          max(abs(p[0] - MIDLINE_X) for p in stray))))
    else:
        checks.append(("PASS", "Beside the ganglion, through nothing",
                       "All %d sampled points are inside the %s foraminal tissue, "
                       "none are in the ganglion and none are in the thecal sac, "
                       "and the lead surface comes within about %.2f mm of the "
                       "ganglion." % (len(pts), target, gap_drg)))

    # Grazing a root sheath is reported but is NOT a failure, and the reason is
    # RADO's own lead: measured against these same meshes it clips three root
    # bodies -- 7 points of a branch root's core, 8 and 56 of two "Middle"
    # layers -- while sitting 100% inside the foraminal tissue and 0% inside the
    # ganglion. RADO's compartments here are not carved out of each other (the
    # foraminal block CONTAINS the root and the ganglion), so a small overlap
    # with a root sleeve is a property of the model, not of the placement. It
    # still matters for meshing, so it is counted and named.
    #
    # Note also that tissue_map.yaml flags the roots' "Middle" layer as possibly
    # CSF rather than nerve; if that is settled the other way, an overlap there
    # means even less than it does now.
    if pierced:
        checks.append(("WARN", "Grazes a nerve-root body",
                       "%d sampled point%s of the lead surface fall inside %d "
                       "nerve-root %s near this foramen (%s). RADO's own lead does "
                       "the same, because the foraminal block is not carved out "
                       "around the roots -- but the mesher will have to resolve the "
                       "overlap, so it is worth seeing before a solve."
                       % (sum(n for _f, n in pierced),
                          "" if sum(n for _f, n in pierced) == 1 else "s",
                          len(pierced), "body" if len(pierced) == 1 else "bodies",
                          ", ".join(fn.replace("T8-10 - ", "")[:38]
                                    for fn, _n in pierced[:3]))))

    ok = not any(level == "FAIL" for level, _t, _m in checks)
    return {"ok": ok, "checks": checks, "numbers": numbers, "points": pts,
            "shapes": shapes, "segments": segments, "params": params}


# --------------------------------------------------------------------------
# validation -- the set, i.e. the checks no single lead can make
# --------------------------------------------------------------------------
def lead_solid(shapes):
    """One lead's segments fused into a single shape, for pairwise distance.

    Fused rather than left as a compound because distToShape on a compound of
    seventeen cylinders does seventeen-squared pair tests, and because a fused
    solid is what common() needs to give a meaningful overlap volume.
    """
    solid = shapes[0][2]
    for _kind, _idx, shape in shapes[1:]:
        solid = solid.fuse(shape)
    return solid


def pair_clearance(a_result, b_result):
    """How close two built leads come, and whether they overlap. A dict.

    The gap is the exact shape-to-shape distance -- not a sampled one -- because
    "do these two 1.3 mm cylinders touch" is a question OCC answers exactly and
    a vertex sample answers only probably. When the gap is zero the overlap is
    quantified two ways: the shared volume, and the fraction of one lead's
    sampled surface that is inside the other, which is the number that reads
    clinically ("a third of lead 2 is inside lead 1").
    """
    solid_a = lead_solid(a_result["shapes"])
    solid_b = lead_solid(b_result["shapes"])
    gap = solid_a.distToShape(solid_b)[0]
    out = {"gap_mm": gap, "overlap_mm3": 0.0, "pct_a_in_b": 0.0}
    if gap > 1e-7:
        return out
    try:
        out["overlap_mm3"] = solid_a.common(solid_b).Volume
    except Exception:
        out["overlap_mm3"] = float("nan")
    pts = a_result.get("points") or sample_points(a_result["shapes"])
    _fc, _Mesh, Part, _msl, _mc, _mf = _geometry_modules()
    n_in = 0
    for p in pts:
        try:
            if solid_b.isInside(Part.Vertex(*p).Point, 1e-6, True):
                n_in += 1
        except Exception:
            break
    out["pct_a_in_b"] = 100.0 * n_in / len(pts) if pts else 0.0
    return out


def validate_set(leads, meshes=None, corridor=None, foramen=None,
                 max_samples=6000, stl_dir=DEFAULT_STL):
    """Validate every lead, then everything only the set can see.

    Returns {"ok", "leads": [result, ...], "pairs": [...], "checks": [...]},
    where "checks" holds the SET-level verdicts -- collisions today -- and each
    lead's own verdicts stay on its own result.
    """
    if meshes is None:
        meshes = load_meshes(stl_dir)
    if corridor is None:
        corridor = json.load(open(DEFAULT_CORRIDOR))
    if foramen is None and any(l["type"] == "drg" for l in leads):
        foramen = load_foramen_corridors()

    results = [validate_lead(lead, meshes=meshes, corridor=corridor,
                             foramen=foramen, max_samples=max_samples,
                             stl_dir=stl_dir) for lead in leads]

    checks, pairs = [], []
    built = [(i, r) for i, r in enumerate(results) if r.get("shapes")]
    if len(built) < 2:
        if len(leads) > 1:
            checks.append(("WARN", "Leads could not be compared",
                           "Fewer than two of the %d leads could be built, so they "
                           "have not been checked against each other." % len(leads)))
    else:
        for ai in range(len(built)):
            for bi in range(ai + 1, len(built)):
                i, a = built[ai]
                j, b = built[bi]
                info = pair_clearance(a, b)
                info.update(a=i, b=j,
                            a_label=leads[i]["label"], b_label=leads[j]["label"])
                pairs.append(info)

        worst = min(pairs, key=lambda p: p["gap_mm"])
        touching = [p for p in pairs if p["gap_mm"] <= 1e-7]
        tight = [p for p in pairs if 1e-7 < p["gap_mm"] < TIGHT_LEAD_GAP_MM]
        if touching:
            first = touching[0]
            checks.append((
                "FAIL", "Two leads occupy the same space",
                "%s and %s intersect: they share %.3f mm3 of volume, and %.1f%% of "
                "the first lead's sampled surface is inside the second. Two leads "
                "cannot both be there. Move them apart -- for a dorsal pair, "
                "increase the difference between their x_offset values until the "
                "gap between them is comfortably more than zero."
                % (first["a_label"], first["b_label"], first["overlap_mm3"],
                   first["pct_a_in_b"])))
        elif tight:
            first = min(tight, key=lambda p: p["gap_mm"])
            checks.append((
                "WARN", "Two leads are very close",
                "%s and %s come within %.2f mm of each other. They do not touch, but "
                "a mesher given a gap that small between two curved conductors will "
                "produce sliver elements there, and the field between them will be "
                "dominated by the gap rather than by the anatomy."
                % (first["a_label"], first["b_label"], first["gap_mm"])))
        else:
            checks.append((
                "PASS", "The leads clear each other",
                "The closest approach between any two of the %d leads is %.2f mm "
                "(%s and %s). None of them intersect."
                % (len(leads), worst["gap_mm"], worst["a_label"], worst["b_label"])))

    ok = all(r["ok"] for r in results) \
        and not any(level == "FAIL" for level, _t, _m in checks)
    return {"ok": ok, "leads": results, "pairs": pairs, "checks": checks}


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------
def report_lines(params, result):
    """One lead's validation as text, for the report file and the macro's panel."""
    n = result["numbers"]
    lines = ["%s -- %s" % (params["label"], "PASS" if result["ok"] else "FAIL"),
             describe(params), ""]
    if params.get("note"):
        lines += [params["note"], ""]
    for level, title, message in result["checks"]:
        lines.append("[%-4s] %s" % (level, title))
        for chunk in wrap(message, 74):
            lines.append("        " + chunk)
    lines.append("")
    lines.append("measurements")
    if params["type"] == "drg":
        if "lead_u" in n:
            lines.append("   lead        %.1f .. %.1f mm from the midline (%.1f mm), "
                         "contacts %.1f .. %.1f"
                         % (n["lead_u"][0], n["lead_u"][1], n["length_mm"],
                            n["contacts_u"][0], n["contacts_u"][1]))
            lines.append("   target      ganglion %s (%s side), spanning %.1f .. %.1f mm "
                         "from the midline at z = %.1f"
                         % (n["target"], n["side"], n["ganglion_u"][0],
                            n["ganglion_u"][1], n["ganglion_z"]))
        if "wall_room_mm" in n:
            lines.append("   foramen     %.2f mm from the lead axis to the wall at "
                         "%.1f mm out; lead %.2f mm; clearance %+.2f mm"
                         % (n["wall_room_mm"], n["wall_room_at_u"],
                            params["diameter"], n["clearance_mm"]))
            lines.append("   containment %.1f%% of %d sampled vertices in foraminal "
                         "tissue, %.1f%% in the ganglion, %.1f%% in dura"
                         % (n["pct_in_foraminal_tissue"], n["sampled_vertices"],
                            n["pct_in_ganglion"], n["pct_in_dura"]))
            lines.append("   nearest     about %.2f mm from the ganglion surface "
                         "(asked for %.2f mm; the lead rides %.2f mm rostral of "
                         "the corridor centre to get it)"
                         % (n["gap_to_ganglion_mm"],
                            n.get("ganglion_clearance_asked_mm", 0.0),
                            n.get("ganglion_standoff_mm", 0.0)))
    else:
        lines.append("   lead        z %.1f .. %.1f mm (%.1f mm), contacts z %.1f .. %.1f"
                     % (n["lead_z"][0], n["lead_z"][1], n["length_mm"],
                        n["contacts_z"][0], n["contacts_z"][1]))
        lines.append("   position    x = %.2f mm (midline %.2f%+.2f)"
                     % (n["x"], MIDLINE_X, float(params["x_offset"])))
        if "channel_min_mm" in n:
            lines.append("   channel     min %.2f mm at z %.1f, median %.2f mm; lead %.2f mm; "
                         "clearance %+.2f mm"
                         % (n["channel_min_mm"], n["channel_min_at_z"],
                            n["channel_median_mm"], params["diameter"],
                            n["clearance_mm"]))
            lines.append("   containment %.1f%% of %d sampled vertices in epidural fat, "
                         "%.1f%% in dura"
                         % (n["pct_in_fat"], n["sampled_vertices"], n["pct_in_dura"]))
    if "centreline_source" in n:
        lines.append("   centreline  %s" % n["centreline_source"])
        lines.append("   local measurement vs stored fit: max %.3f mm apart"
                     % n["centreline_vs_stored_max_mm"])
    return lines


def report_set_lines(name, leads, outcome):
    """The whole configuration's validation as text."""
    lines = ["%s -- %s" % (name, "PASS" if outcome["ok"] else "FAIL"),
             describe_set(leads), ""]
    for lead, result in zip(leads, outcome["leads"]):
        lines.append("-" * 72)
        lines += report_lines(lead, result)
        lines.append("")
    if outcome["checks"]:
        lines.append("-" * 72)
        lines.append("the set as a whole")
        for level, title, message in outcome["checks"]:
            lines.append("[%-4s] %s" % (level, title))
            for chunk in wrap(message, 74):
                lines.append("        " + chunk)
    if outcome["pairs"]:
        lines.append("")
        lines.append("closest approach between leads")
        for p in sorted(outcome["pairs"], key=lambda q: q["gap_mm"]):
            lines.append("   %-22s %-22s %7.2f mm%s"
                         % (p["a_label"], p["b_label"], p["gap_mm"],
                            "   INTERSECTING" if p["gap_mm"] <= 1e-7 else ""))
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
# validating the DRG trajectory against RADO's own lead
# --------------------------------------------------------------------------
def compare_rado(params, result, stl_dir=DEFAULT_STL):
    """Measure a generated DRG lead against RADO's own 4-contact lead. Lines.

    This is the only independent check the DRG trajectory has. RADO placed that
    lead by hand in SolidWorks, on the left level-3 ganglion, and it is not
    derived from anything in this pipeline -- so if a lead generated from the
    ray-cast foraminal corridor lands where RADO's is, the corridor is right.
    Contact centres are compared pairwise in the order they run out along the
    root, because RADO's four are NOT numbered in spatial order (by x they run
    3, 1, 2, 4) and comparing them by index would be meaningless.
    """
    _fc, Mesh, _p, _msl, _mc, mf = _geometry_modules()
    lines = []
    if params["type"] != "drg":
        return ["--compare-rado only means anything for a DRG lead."]
    if params["target"] != RADO_DRG_TARGET:
        lines.append("NOTE: RADO's own lead sits on ganglion %s, and this one targets "
                     "%s, so the comparison below is between different foramina and "
                     "says nothing about the trajectory."
                     % (RADO_DRG_TARGET, params["target"]))

    rado = []
    for fn in RADO_DRG_STLS:
        mesh = Mesh.Mesh(os.path.join(stl_dir, fn))
        cx, cy, cz = mf.area_centroid(mesh)
        rado.append((abs(cx - MIDLINE_X), cy, cz, fn))
    rado.sort()

    mine = []
    for (kind, idx, shape), (_k, _i, ts, ln) in zip(result["shapes"],
                                                    result["segments"]):
        if kind != "contact":
            continue
        c = shape.CenterOfMass
        mine.append((abs(c.x - MIDLINE_X), c.y, c.z, "generated contact %d" % idx))
    mine.sort()

    # -- the trajectory question, which is separate from the placement one --
    #
    # Where the contacts sit along the root is a PARAMETER (lateral_offset), so
    # comparing contact centres mixes "is the curve right" with "did we choose
    # the same spot on it". The curve on its own is compared here: at each u
    # where RADO put a contact, how far is the generated corridor's own point
    # from RADO's contact centre, across the corridor rather than along it?
    corridor = None
    try:
        corridor = corridor_for(params, load_foramen_corridors())
    except ValueError:
        pass
    if corridor is not None:
        cy = corridor["y_poly_coef_high_to_low"]
        cz = corridor["z_poly_coef_high_to_low"]
        standoff = float((result.get("numbers") or {}).get("ganglion_standoff_mm", 0.0))
        lines += ["the TRAJECTORY, independent of where the array was placed on it",
                  "(the fitted foraminal corridor, plus this lead's derived %.2f mm "
                  "rostral standoff, sampled at RADO's own contact positions)"
                  % standoff, "",
                  "   %-8s %-16s %-16s %s"
                  % ("at u", "RADO contact", "corridor here", "across-corridor gap"),
                  "   %-8s %7s %7s  %7s %7s  %7s %7s %7s"
                  % ("", "y", "z", "y", "z", "dy", "dz", "total")]
        across = []
        for u, ry, rz, _fn in rado:
            gy = mf.polyval(cy, u) + float(params.get("y_offset", 0.0))
            gz = mf.polyval(cz, u) + float(params.get("z_offset", 0.0)) + standoff
            d = ((gy - ry) ** 2 + (gz - rz) ** 2) ** 0.5
            across.append((u, d))
            lines.append("   %-8.2f %7.2f %7.2f  %7.2f %7.2f  %+7.2f %+7.2f %7.2f"
                         % (u, ry, rz, gy, gz, gy - ry, gz - rz, d))
        over = [d for u, d in across if u >= corridor["ganglion_u"][0]]
        lines += ["",
                  "   across-corridor gap: mean %.2f mm over all four, %.2f mm over "
                  "the %d that lie at or beyond the ganglion"
                  % (sum(d for _u, d in across) / len(across),
                     (sum(over) / len(over)) if over else float("nan"), len(over)),
                  ""]

    lines += ["the CONTACTS, i.e. this configuration's actual placement",
              "(paired in order along the root, medial first; RADO's own "
              "numbering is not spatial)", "",
              "   %-8s %-24s %-24s %s" % ("", "RADO", "generated", "difference"),
              "   %-8s %7s %7s %7s  %7s %7s %7s  %6s %6s %6s"
              % ("contact", "u", "y", "z", "u", "y", "z", "du", "dy", "dz")]
    n = min(len(rado), len(mine))
    dists = []
    for i in range(n):
        r, m = rado[i], mine[i]
        d = ((r[0] - m[0]) ** 2 + (r[1] - m[1]) ** 2 + (r[2] - m[2]) ** 2) ** 0.5
        dists.append(d)
        lines.append("   %-8d %7.2f %7.2f %7.2f  %7.2f %7.2f %7.2f  %+6.2f %+6.2f %+6.2f"
                     % (i + 1, r[0], r[1], r[2], m[0], m[1], m[2],
                        m[0] - r[0], m[1] - r[1], m[2] - r[2]))
    if len(rado) != len(mine):
        lines.append("   (RADO has %d contacts, this lead has %d; the first %d are "
                     "compared)" % (len(rado), len(mine), n))
    if dists:
        lines += ["",
                  "   contact-centre separation: mean %.2f mm, worst %.2f mm"
                  % (sum(dists) / len(dists), max(dists))]
    return lines


# --------------------------------------------------------------------------
# outputs: STLs on disk, and a disposable preview in a live document
# --------------------------------------------------------------------------
def lead_dir_name(params):
    """The subdirectory one lead of a multi-lead set is exported into."""
    if params["type"] == "drg":
        return "lead%d_drg_%s" % (params["index"], params["target"])
    return "lead%d_%s" % (params["index"], params["type"])


def export_set(name, leads, results, out_root=DEFAULT_OUT, tag=None):
    """Write the STLs for a whole configuration. Returns (dir, [relative files]).

    Layout, and it is deliberately not uniform:

        one lead   generated_leads/<tag>/SCS Lead Electrode 1.stl, ...
        several    generated_leads/<tag>/lead1_dorsal/SCS Lead Electrode 1.stl, ...

    A single-lead configuration keeps the flat layout it has always had, because
    dorsal_z110_8c and ventral_z110_8c are committed at those exact paths and
    are the regression test for this whole module; moving them would make
    "reproduces its own STLs byte for byte" untestable. A set gets one
    subdirectory per lead, plus a leads.txt manifest, because the FILENAMES
    cannot carry the distinction: they are RADO's own -- "SCS Lead Electrode
    1.stl", "SCS Lead Insulator.stl" -- and ../ansys/tissue_map.yaml matches the
    object Names those sanitise to, so renaming them costs the contacts their
    silver and the insulator its material. The directory carries it instead.

    Importing two leads into one FreeCAD document therefore gives the second
    lead's bodies uniquified Names (SCS_Lead_Electrode_001..; FreeCAD strips the
    trailing digits before appending its counter). tissue_map.yaml still matches
    them -- it keys on the "SCS Lead Electrode" prefix -- but which contact
    belongs to which lead is then only recoverable from leads.txt, so import one
    lead at a time and relabel, or read the manifest.

    Nothing here opens or saves a project document; the scratch document exists
    only because Mesh.export() takes document objects.
    """
    FreeCAD, Mesh, _p, msl, _mc, _mf = _geometry_modules()
    root = os.path.join(out_root, tag or name)
    if not os.path.isdir(root):
        os.makedirs(root)
    several = len(leads) > 1

    doc = FreeCAD.newDocument("scs_lead_export")
    written = []
    manifest = ["%s -- %s" % (tag or name, describe_set(leads)), ""]
    try:
        for params, result in zip(leads, results):
            outdir = os.path.join(root, lead_dir_name(params)) if several else root
            if not os.path.isdir(outdir):
                os.makedirs(outdir)
            manifest += ["%s  ->  %s" % (params["label"],
                                         os.path.relpath(outdir, root)),
                         "   %s" % describe(params)]
            for kind, idx, shape in result["shapes"]:
                if kind != "contact":
                    continue
                obj = doc.addObject("Part::Feature", "contact%d" % idx)
                obj.Shape = shape
                path = os.path.join(outdir, "SCS Lead Electrode %d.stl" % idx)
                Mesh.export([obj], path)
                written.append(os.path.relpath(path, root))
            ins = doc.addObject("Part::Feature", "insulator")
            ins.Shape = msl.fuse_insulator(result["shapes"])
            path = os.path.join(outdir, "SCS Lead Insulator.stl")
            Mesh.export([ins], path)
            written.append(os.path.relpath(path, root))
            manifest.append("")
    finally:
        FreeCAD.closeDocument(doc.Name)
    if several:
        with open(os.path.join(root, "leads.txt"), "w", encoding="utf-8") as fh:
            fh.write("\n".join(manifest) + "\n")
    return root, written


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


def preview_set(doc, name, leads, results, colours=None):
    """Put a whole configuration into a live document as a disposable preview.

    Replaces any previous preview, so exploring configurations does not silt the
    tree up. The bodies are Part::Feature -- parametric shapes, a few kB -- not
    imported meshes, so a preview costs nothing and leaves nothing on disk. Which
    is the whole point: the configuration is the artefact, the document is not.

    NAMES ARE PREFIXED PER LEAD, and that is not cosmetic. FreeCAD uniquifies a
    colliding Name by STRIPPING ITS TRAILING DIGITS and appending a 3-digit
    counter, so a second lead's "SCS_Preview_Contact1" would come back as
    "SCS_Preview_Contact001" and the two leads' contacts would alphabetise into
    each other -- exactly the confusion that once made one 8-contact lead read as
    a 12-contact one in the tree. Giving each lead its own
    "SCS_Preview_LeadN_Contact_MM" prefix means no Name ever collides, so none
    is ever rewritten, and the lead number and contact number in the Name both
    mean what they say. Contact numbers are zero-padded so they sort in physical
    order (contact 01 is the caudal end of an epidural lead, and the MEDIAL end
    of a DRG lead -- the end nearest the cord).

    THIS FUNCTION DOES NOT SAVE, and must not learn to. It refuses outright on a
    document that has a FileName when there is no GUI, because that is the shape
    of the accident that destroys the model: a headless run leaves the real
    document dirty, someone saves it, and every ShapeAppearance in it is gone.
    See make_tissue_groups.py.
    """
    FreeCAD, _Mesh, _p, msl, _mc, _mf = _geometry_modules()
    if not FreeCAD.GuiUp and getattr(doc, "FileName", ""):
        raise RuntimeError(
            "refusing to preview into %s with no GUI running: a headless save "
            "would strip the document's colours. Preview from a running FreeCAD, "
            "or export STLs instead." % os.path.basename(doc.FileName))

    clear_preview(doc)
    group = doc.addObject("App::DocumentObjectGroup", PREVIEW_GROUP)
    group.Label = "PREVIEW: %s" % name

    made = []
    for params, result in zip(leads, results):
        if not result.get("shapes"):
            continue
        n = params["index"]
        for kind, idx, shape in result["shapes"]:
            if kind != "contact":
                continue
            obj = doc.addObject("Part::Feature",
                                "SCS_Preview_Lead%d_Contact_%02d" % (n, idx))
            obj.Shape = shape
            obj.Label = "%s — contact %02d" % (params["label"], idx)
            made.append(("electrode_contact", obj))
        ins = doc.addObject("Part::Feature", "SCS_Preview_Lead%d_Insulator" % n)
        ins.Shape = msl.fuse_insulator(result["shapes"])
        ins.Label = "%s — insulator" % params["label"]
        made.append(("lead_insulation", ins))

    group.Group = [obj for _t, obj in made]
    if colours:
        for tissue, obj in made:
            rgb = colours.get(tissue)
            if rgb and obj.ViewObject is not None:
                paint(obj.ViewObject, rgb)
    doc.recompute()
    return group, [obj for _t, obj in made]


RADO_REFERENCE_GROUP = "SCS_RadoLeadReference"
RADO_LEAD_STLS = RADO_DRG_STLS + ["SCS Lead Insulator.stl"]


def preview_rado_lead(doc, stl_dir=DEFAULT_STL, colours=None):
    """Drop RADO's own 4-contact DRG lead into the document, from its own STLs.

    This is what makes removing those five bodies from NBF_RADO-SCS.FCStd a
    change of storage rather than a loss. The five files are RADO's, unmodified,
    tracked in STL_files/, and already in the model's coordinate frame -- they
    are imported with NO transform, exactly as the published model has them --
    so the geometry is one button away whenever it is wanted, and is NOT sitting
    in the anatomy document acting as a second conductor in every field solve.

    Its own group, not the lead preview's, so previewing a generated lead beside
    RADO's does not delete RADO's -- looking at the two together is the whole
    point. Cleared with clear_rado_lead().
    """
    FreeCAD, Mesh, _p, _msl, _mc, _mf = _geometry_modules()
    if not FreeCAD.GuiUp and getattr(doc, "FileName", ""):
        raise RuntimeError(
            "refusing to import into %s with no GUI running: a headless save "
            "would strip the document's colours."
            % os.path.basename(doc.FileName))
    clear_rado_lead(doc)
    before = len(doc.Objects)
    for fn in RADO_LEAD_STLS:
        Mesh.insert(os.path.join(stl_dir, fn), doc.Name)
    added = doc.Objects[before:]
    group = doc.addObject("App::DocumentObjectGroup", RADO_REFERENCE_GROUP)
    group.Label = "REFERENCE: RADO's own 4-contact DRG lead"
    for fn, obj in zip(RADO_LEAD_STLS, added):
        obj.Label = ("RADO DRG lead insulator" if "Insulator" in fn
                     else "RADO DRG lead contact %s" % os.path.splitext(fn)[0][-1])
        if colours and obj.ViewObject is not None:
            rgb = colours.get("lead_insulation" if "Insulator" in fn
                              else "electrode_contact")
            if rgb:
                paint(obj.ViewObject, rgb)
    group.Group = list(added)
    doc.recompute()
    return group, list(added)


def clear_rado_lead(doc):
    """Remove the RADO reference lead overlay. Returns how many objects went."""
    group = doc.getObject(RADO_REFERENCE_GROUP)
    if group is None:
        return 0
    names = [o.Name for o in group.Group]
    for name in names:
        doc.removeObject(name)
    doc.removeObject(group.Name)
    doc.recompute()
    return len(names) + 1


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
    ap.add_argument("--foramen", default=DEFAULT_FORAMEN)
    ap.add_argument("--stl-dir", default=DEFAULT_STL)
    ap.add_argument("--out", default=DEFAULT_OUT)
    ap.add_argument("--report", default=DEFAULT_REPORT)
    ap.add_argument("--name", help="a configuration from lead_configs.yaml")
    ap.add_argument("--all", action="store_true", help="every configuration in the file")
    ap.add_argument("--list", action="store_true", help="list them; no FreeCAD needed")
    ap.add_argument("--show", action="store_true",
                    help="resolve to numbers and stop; no FreeCAD needed")
    ap.add_argument("--validate", action="store_true", help="measure the fit")
    ap.add_argument("--compare-rado", action="store_true",
                    help="measure a DRG lead against RADO's own 4-contact lead")
    ap.add_argument("--export", action="store_true",
                    help="validate, then write STLs if it passes")
    ap.add_argument("--force", action="store_true",
                    help="export even if validation failed (it will not fit)")
    ap.add_argument("--tag", help="output directory name; defaults to the config name")
    ap.add_argument("--max-samples", type=int, default=6000)
    for flag, kind in (("--type", str), ("--side", str), ("--contacts", int),
                       ("--contact-length", float), ("--gap", float),
                       ("--diameter", float), ("--x-offset", float),
                       ("--z-center", float), ("--tail", float), ("--level", str),
                       ("--target", str), ("--lateral-offset", float),
                       ("--ganglion-clearance", float),
                       ("--y-offset", float), ("--z-offset", float)):
        ap.add_argument(flag, type=kind, default=None, help="override the configuration")
    args = ap.parse_args(script_args())

    out = []

    def say(fmt, *a):
        out.append(fmt % a if a else fmt)

    cfg = parse_config_file(args.configs)
    overrides = {k: getattr(args, k.replace("-", "_"))
                 for k in ("type", "side", "contacts", "contact_length", "gap",
                           "diameter", "x_offset", "z_center", "tail", "level",
                           "target", "lateral_offset", "ganglion_clearance",
                           "y_offset", "z_offset")}

    if args.list:
        say("configurations in %s", os.path.relpath(args.configs, REPO))
        say("")
        for name in config_names(cfg):
            try:
                leads = resolve_set(name, cfg)
            except ValueError as exc:
                say("%-24s CANNOT RESOLVE: %s", name, exc)
                continue
            say("%-24s %s", name, describe_set(leads))
            if len(leads) > 1:
                for lead in leads:
                    say("%-24s    %s", "", describe(lead))
        say("")
        say("at most %s leads of any one type per configuration",
            cfg.get("max_leads_per_type", DEFAULT_MAX_LEADS_PER_TYPE))
        say("")
        say("vertebral levels (z of the body centre, mm)")
        for level, z in cfg["levels"].items():
            say("   %-14s %8.2f", level, float(z))
        if os.path.exists(args.foramen):
            foramen = json.load(open(args.foramen))
            say("")
            say("ganglia a DRG lead can target (rostral to caudal within each side)")
            for tag in sorted(foramen["ganglia"]):
                g = foramen["ganglia"][tag]
                say("   %-6s %-6s side, core at z %7.2f, foramen measurable "
                    "%.1f-%.1f mm from the midline",
                    tag, g["side"], g["ganglion_centroid"][2],
                    g["usable_u"][0], g["usable_u"][1])
        else:
            say("")
            say("no foraminal corridor measurements yet -- run "
                "`freecadcmd src/freecad/measure_foramen.py` before using a DRG lead")
        say("")
        say("RESULT: %d configurations", len(cfg["configs"]))
        write_report(args.report, out)
        return 0

    names = config_names(cfg) if args.all else [args.name]
    if names == [None]:
        say("Nothing to do: give --name NAME, or --all, or --list.")
        write_report(args.report, out)
        return 2

    # Fail here rather than three functions deep with an ImportError traceback:
    # measuring anything needs FreeCAD, and the fix is a different interpreter,
    # not a different configuration.
    if args.validate or args.export or args.compare_rado:
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

    meshes = load_meshes(args.stl_dir) if (args.validate or args.export
                                           or args.compare_rado) else None
    corridor = json.load(open(args.corridor)) if meshes is not None else None
    foramen = (load_foramen_corridors(args.foramen)
               if meshes is not None and os.path.exists(args.foramen) else None)

    rc = 0
    for name in names:
        try:
            leads = resolve_set(name, cfg, overrides)
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
        say("%s", describe_set(leads))
        if args.show:
            for i, lead in enumerate(leads, start=1):
                say("   -- %s", lead["label"])
                for key in sorted(PARAM_TYPES):
                    if key in lead:
                        say("      %-16s %s", key, lead[key])
                if "x" in lead:
                    say("      %-16s %.2f", "x (absolute)", lead["x"])
            say("")
            continue

        try:
            outcome = validate_set(leads, meshes=meshes, corridor=corridor,
                                   foramen=foramen, max_samples=args.max_samples,
                                   stl_dir=args.stl_dir)
        except ValueError as exc:
            say("   could not validate: %s", exc)
            say("")
            rc = 2
            continue
        say("")
        for line in report_set_lines(name or "custom", leads, outcome)[2:]:
            say("%s", line)
        if not outcome["ok"]:
            rc = 1

        if args.compare_rado:
            for lead, result in zip(leads, outcome["leads"]):
                if lead["type"] != "drg" or not result.get("shapes"):
                    continue
                say("")
                say("-" * 72)
                for line in compare_rado(lead, result, args.stl_dir):
                    say("%s", line)

        if args.export:
            if not outcome["ok"] and not args.force:
                say("")
                say("NOT EXPORTED: this configuration does not fit. Fix it, or pass "
                    "--force if you know what you are doing.")
            else:
                outdir, files = export_set(name or "custom", leads, outcome["leads"],
                                           args.out, args.tag)
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
