#!/usr/bin/env python3
"""Build and VALIDATE one SCS lead from its parameters. No named configurations.

WHAT A CLINICIAN ASKS FOR, AND WHAT THIS TURNS IT INTO
------------------------------------------------------
The clinical description of a lead is "8 contacts, 3 mm long, 1 mm apart,
dorsal, at T10, 1 mm left of midline", or "four contacts on the left third
ganglion". None of that is geometry. This module is the layer that turns it into
geometry and then, more importantly, tells you whether the answer is
anatomically possible in THIS model:

    parameters         ->  resolve_lead()  ->  one lead, a flat set of numbers
    (flags, or the         resolve_leads()     a set of them, capped at
     panel's spinboxes)                        max_leads_per_type of each type
                       ->  validate_set()  ->  PASS/WARN/FAIL per lead AND for
                                               the set as a whole
                       ->  export_set()    ->  generated_leads/<name>/*.stl
                       ->  preview_set()   ->  disposable bodies in the open
                                               document

src/freecad/lead_designer.FCMacro is the same thing with spinboxes on it.

THERE IS NO CATALOGUE OF PLACEMENTS
-----------------------------------
There used to be a `configs:` section in the YAML, one entry per named
placement, and it grew every time a lead was positioned. It is gone. A lead is
its parameters; naming one costs an invented name, a paragraph explaining it,
and a standing obligation to keep the two in step with the geometry. What
survives is lead_defaults.yaml -- the starting values, the type-specific default
blocks, the vertebral levels and max_leads_per_type -- which is small and stable
because it describes what a lead IS, not where anyone once put one.

Two placements ARE pinned, in src/freecad/test_regression.py rather than in the
YAML: the dorsal and ventral 8-contact leads whose 18 STLs are committed under
generated_leads/. That test holds their parameters explicitly and asserts the
exported files come back byte for byte, which is what has twice caught this
module drifting. It is a test, not a catalogue: nothing builds from it.

THREE LEAD TYPES, AND TWO OF THEM ARE BUILT WHILE THE THIRD IS COPIED
---------------------------------------------------------------------
    dorsal / ventral    BUILT. Swept along the epidural canal centreline y(z)
                        measured at the requested x by measure_corridor.py. The
                        lead runs rostro-caudally; z is the sweep parameter and
                        x is fixed. Fully parametric: contact count, contact
                        length, gap, diameter and tail are all yours.

    drg                 COPIED. RADO's own four-contact lead -- the five STLs in
                        STL_files/, hand-placed in SolidWorks against the left
                        third ganglion -- transformed rigidly onto whichever
                        ganglion is named. See rado_drg_lead.py for the frame
                        the transform is built on, and for why going to the
                        right side is a REFLECTION rather than a rotation.
                        Nothing about its hardware is parametric: a rigid map
                        cannot resize a contact.

THE DRG LEAD USED TO BE SWEPT TOO, AND THAT WAS WRONG. Sweeping it along the
foraminal corridor produced a lead that fitted the foramen without ever being
the lead RADO drew -- its sheath curvature was whatever the sweep gave, and it
had to be lifted off the corridor centre by a solved-for standoff to stop it
cutting into the ganglion it was supposed to lie beside. RADO's lead needs
neither correction, because RADO placed it correctly to begin with. So the
foraminal corridor stopped being a construction and became what it should have
been all along: a measurement, used to say how much room there is at a target
and to cross-check the frame.

The three are not settings of one thing, so their parameters do not overlap and
resolve_lead() REJECTS rather than ignores one that means nothing for the type.
`z_center` on a DRG lead, or `target` on a dorsal one, is an error, on the same
principle that a typo'd `diamter: 2.0` is an error rather than a silently kept
default -- and so, now, is `contacts` on a DRG lead, which is not a typo but a
request for something this type of lead cannot do.

WHY LATERAL OFFSETS RE-MEASURE THE CENTRELINE
---------------------------------------------
epidural_corridor.json holds the centreline y(z) measured at the MIDLINE. The
channel centre moves as you go lateral -- about 1.1 mm of y at 4 mm off the
midline -- so a lead swept along the midline curve at x = midline + 4 sits about
half its clearance off centre before it has done anything else. Any lead with a
non-zero x_offset therefore gets its own centreline, ray cast and fitted at ITS
x over ITS z span. At x_offset = 0 this reproduces the stored polynomial; the
report prints the agreement so the two can never silently diverge.

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

A DRG lead asks different questions, for two reasons. The foramen is not a
void -- RADO models the tissue filling each neuroforamen as one SOLID block, and
a percutaneous DRG lead is threaded through foraminal fat rather than dropped
into a gap. And the lead is copied rather than built, so two of the questions a
swept lead has to answer do not arise: whether the contacts came out the right
size, and whether the trajectory matches RADO's. They are RADO's. So:

    placement   the transform is rigid, its determinant is +1 within a side and
                -1 across sides, and every transformed body came back watertight
                with OUTWARD normals -- a cross-side copy is a reflection and
                inverts them unless the winding is reversed with it
    clearance   every sampled vertex OUTSIDE the target ganglion and OUTSIDE the
                thecal sac. Either one is a FAIL naming which, because a lead
                that passes through the ganglion it is meant to stimulate is not
                a placement, it is an injury
    containment every sampled vertex inside the foraminal block -- and where one
                is not, WHAT it is in. Out of the DISTAL MOUTH into space no
                compartment occupies is what a real lead's tip does, and warns;
                into a pedicle or a disc fails
    standoff    how far the lead surface sits off the ganglion surface, against
                the same measurement for RADO's own lead at its own ganglion.
                This is the number that says the copy still "interfaces smoothly
                with the corresponding DRG"
    foramen     how much channel there is at this target, from measure_foramen's
                stations. It no longer decides anything; it is what EXPLAINS one
                target being tighter than another

AND ONE CHECK THE SET OWNS RATHER THAN ANY LEAD
-----------------------------------------------
Two leads can each pass every check above and still be un-implantable, because
nothing in a per-lead check looks at the OTHER lead. Two 1.30 mm dorsal leads
1.0 mm apart both fit the channel; they simply cannot both be there. So
validate_set() measures every pair -- for two swept leads the exact
shape-to-shape distance and, when that is zero, the volume they share and how
much of one lead's sampled surface is inside the other; for any pair involving a
copied DRG lead the nearest point-to-triangle approach between the two surfaces,
because there is no solid to hand OCC. pair_clearance() records which it did.

A lead that fails is REPORTED, in clinical terms, and not exported. It is not
carved into the anatomy to make it fit.

ONE MODEL DOCUMENT, AND IT HOLDS NO LEAD
----------------------------------------
NBF_RADO-SCS.FCStd holds anatomy and nothing else -- 240 bodies in 11 groups,
RADO's own 4-contact DRG lead removed by drop_rado_lead.py. There is no .FCStd
per placement: the durable artefact is generated_leads/<name>/, and the bodies
in a FreeCAD document are a DISPOSABLE PREVIEW:

    preview_set(doc, leads, results)  adds every lead to the open document
                                      inside one group labelled "12 SCS Leads",
                                      replacing any previous preview
    clear_preview(doc)                deletes that group and everything in it

Try ten arrangements, export the one you want, and close the document WITHOUT
SAVING.

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
Resolving parameters to numbers needs no FreeCAD at all:

    python3 src/freecad/build_lead_config.py --level T11 --x-offset 1.5 --show

Measuring and building does:

    BUILD_LEAD_CONFIG_ARGS="--level T11 --x-offset 1.5 --validate" \\
        freecadcmd src/freecad/build_lead_config.py

    BUILD_LEAD_CONFIG_ARGS="--type drg --target L3 --validate" \\
        freecadcmd src/freecad/build_lead_config.py

    BUILD_LEAD_CONFIG_ARGS="--type ventral --z-center 110.432 --export \\
        --tag ventral_T11" freecadcmd src/freecad/build_lead_config.py

    # RADO's own lead, on its own ganglion, checked against the five STLs it
    # was copied from. The transform is the identity here, so every body must
    # come back exactly where RADO left it -- the DRG path's strongest check:
    BUILD_LEAD_CONFIG_ARGS="--type drg --target L3 --compare-rado" \\
        freecadcmd src/freecad/build_lead_config.py

    # the same lead on the right side: a reflection, not a rotation
    BUILD_LEAD_CONFIG_ARGS="--type drg --target R2 --validate" \\
        freecadcmd src/freecad/build_lead_config.py

Anything not given falls back to lead_defaults.yaml. One lead per run on the
command line; the panel is what builds a SET of them, through the same
resolve_leads() and validate_set().

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
DEFAULTS_YAML = os.path.join(HERE, "lead_defaults.yaml")
DEFAULT_CORRIDOR = os.path.join(HERE, "epidural_corridor.json")
DEFAULT_FORAMEN = os.path.join(HERE, "foraminal_corridors.json")
DEFAULT_STL = os.path.join(REPO, "STL_files")
DEFAULT_OUT = os.path.join(HERE, "generated_leads")
DEFAULT_REPORT = os.path.join(HERE, "build_lead_config_report.txt")

MIDLINE_X = 56.60
EPIDURAL_STL = "T8-10 - neuro_EpiduralSpace-1.STL"
DURA_STL = "T8-10 - neuro_Meninges-1.STL"

# RADO's own 4-contact DRG lead, as it ships in STL_files/. Not part of any
# generated lead -- it is the independent reference the DRG trajectory is
# checked against by compare_rado(); the parameters that reproduce it are in
# this module's docstring, under USAGE.
RADO_DRG_STLS = ["SCS Lead Electrode %d.stl" % i for i in (1, 2, 3, 4)]
RADO_DRG_TARGET = "L3"

# The group every previewed body lives in -- ALL of them, however many leads,
# contacts and insulators there are. Matched on Name, not Label, so retitling it
# in the tree does not orphan the preview.

PREVIEW_GROUP = "SCS_LeadPreview"
PREVIEW_GROUP_LABEL = "12 SCS Leads"

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
# "diamter: 2.0" in lead_defaults.yaml that quietly kept the 1.3 mm default
# would be the worst kind of bug here.
PARAM_TYPES = {
    "type": str,
    "contacts": int,
    "contact_length": float,
    "gap": float,
    "diameter": float,
    "tail": float,
    "note": str,
    "label": str,
    # dorsal / ventral only
    "side": str,
    "x_offset": float,
    "z_center": float,
    "level": str,
    # DRG only
    "target": str,
    "lateral_offset": float,
    "y_offset": float,
    "z_offset": float,
}

# Which parameters each type accepts. The split is the point: a DRG lead has no
# z_center and a dorsal lead has no target, and saying so is what stops a
# configuration from carrying a number that quietly does nothing.
COMMON_PARAMS = {
    "type",
    "contacts",
    "contact_length",
    "gap",
    "diameter",
    "tail",
    "note",
    "label",
}
EPIDURAL_PARAMS = {"side", "x_offset", "z_center", "level"}
DRG_PARAMS = {"target", "lateral_offset", "y_offset", "z_offset"}

# The hardware parameters a DRG lead NO LONGER TAKES. A DRG lead is not swept
# any more -- it is a rigidly transformed copy of the five STLs RADO ships (see
# rado_drg_lead.py), so its contact count, contact length, gap, diameter and
# tail are properties of those meshes and cannot be dialled. They are still
# COMMON_PARAMS, because a dorsal or ventral lead still owns every one of them;
# what changes for a DRG lead is that they are
#
#   - not inherited from lead_defaults.yaml's `defaults:` block, which would
#     otherwise silently attach an 8-contact 1.30 mm spec to a lead that has
#     four 1.25 mm contacts;
#   - REJECTED, not ignored, if anything sets one explicitly, on the same
#     principle as every other parameter that means nothing for its type;
#   - filled back in by MEASURING RADO's meshes, so describe() and the reports
#     can still say what the hardware is without any of it being a claim.
DRG_FIXED_PARAMS = {"contacts", "contact_length", "gap", "diameter", "tail"}

TYPES = ("dorsal", "ventral", "drg")
EPIDURAL_TYPES = ("dorsal", "ventral")

DEFAULT_MAX_LEADS_PER_TYPE = 2


_RADO_HARDWARE = {}


def rado_drg_hardware(stl_dir=DEFAULT_STL):
    """RADO's DRG hardware as a parameter dict, MEASURED off the STLs.

    contacts / contact_length / gap / diameter / tail, filled in for a DRG lead
    in place of the ones it is no longer allowed to set. Measured rather than
    written down, so these cannot drift from the meshes they describe; cached
    per STL directory because resolve_lead() is called once per lead per
    validation and the measurement reads 33000 facets.

    `tail` comes back CLAMPED AT ZERO, and the clamp is worth explaining.
    RADO's end contacts overhang the sheath rather than sitting inboard of it,
    so the measured tails are NEGATIVE (-0.23 mm medially, -0.62 mm laterally).
    There is no such thing as a negative length of plain insulator, and
    check_sane rightly refuses one, so the parameter reads 0.0 -- "no tail" --
    and the true signed overhangs are reported by the fit check, where they are
    a fact about the geometry rather than a number pretending to be a setting.
    """
    if stl_dir not in _RADO_HARDWARE:
        import rado_drg_lead as rdl
        spec = rdl.hardware_spec(rdl.load_rado_bodies(stl_dir))
        _RADO_HARDWARE[stl_dir] = {
            "contacts": int(spec["contacts"]),
            "contact_length": round(float(spec["contact_length"]), 3),
            "gap": round(float(spec["gap"]), 3),
            "diameter": round(float(spec["diameter"]), 3),
            "tail": max(0.0, round(float(spec["tail"]), 3)),
        }
    return dict(_RADO_HARDWARE[stl_dir])


def is_built(result):
    """Did this lead produce geometry? True for either representation.

    A dorsal or ventral lead is a list of swept Part shapes; a DRG lead is a
    list of transformed triangle arrays. Everything that only needs to know
    "is there something to preview, export or collide" asks this instead of
    testing for one of them and quietly skipping the other.
    """
    return bool(result.get("shapes")) or bool(result.get("tris"))


def is_copy(result):
    """True for a DRG lead, which is copied geometry rather than swept."""
    return bool(result.get("tris"))


def is_freecad_interpreter():
    """True under freecadcmd or the GUI console; see apply_colors.py."""
    return os.path.basename(sys.argv[0]).lower().startswith("freecad")


def script_args():
    """Arguments meant for this script, under either interpreter."""
    if is_freecad_interpreter():
        return shlex.split(os.environ.get("BUILD_LEAD_CONFIG_ARGS", ""))
    return sys.argv[1:]


# --------------------------------------------------------------------------
# lead_defaults.yaml
# --------------------------------------------------------------------------
def parse_defaults_file(path=None):
    """Read lead_defaults.yaml.

    Returns the top-level mapping: scalars like `max_leads_per_type` as scalars,
    and `defaults`, `epidural_defaults`, `drg_defaults` and `levels` as mappings
    of scalars. That is the whole shape of the file -- there is no list in it and
    no nesting deeper than two, because it holds STARTING VALUES and not
    placements. The `configs:` catalogue this parser used to also read is gone;
    see the module docstring.

    A hand-rolled parser, for the reason requirements.txt gives: PyYAML is
    deliberately not a dependency of this directory, so the scripts run on a bare
    interpreter with nothing installed. It accepts exactly the subset
    lead_defaults.yaml uses and nothing else -- no flow style, no anchors, no
    multi-line scalars and no lists. Anything it cannot parse RAISES rather than
    being skipped, because a silently dropped parameter is a wrong lead.

    The grammar, by indent:

        0   key: value          a top-level scalar
        0   key:                a section
        2   key: value          a scalar in the section
    """
    path = path or DEFAULTS_YAML
    top, section = {}, None
    for n, raw in enumerate(open(path), 1):
        line = raw.rstrip("\n")
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        indent = len(line) - len(line.lstrip())
        body = line.lstrip()

        m = re.match(r"^([A-Za-z_][A-Za-z0-9_.\-]*)\s*:\s*(.*)$", body)
        if not m:
            raise ValueError("%s:%d: cannot parse %r" % (path, n, line))
        key, rest = m.group(1), strip_comment(m.group(2))

        if indent == 0:
            if rest:
                top[key] = scalar(rest)
                section = None
            else:
                section = {}
                top[key] = section
        elif indent == 2:
            if section is None:
                raise ValueError("%s:%d: %r is outside any section" % (path, n, key))
            if rest == "":
                raise ValueError(
                    "%s:%d: %r has no value -- lead_defaults.yaml "
                    "holds starting values, not nested entries" % (path, n, key)
                )
            section[key] = scalar(rest)
        else:
            raise ValueError("%s:%d: unexpected indent %d" % (path, n, indent))

    for required in ("defaults", "levels"):
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
# resolving parameters into leads
# --------------------------------------------------------------------------
def resolve_leads(cfg, per_lead, overrides=None):
    """Flatten a list of per-lead parameter dicts into resolved leads.

    `per_lead` is one dict per lead -- the panel's tabs, or a single dict of CLI
    flags. `overrides` is applied on top of EVERY lead, which is how one number
    is swept across a pair from the command line.

    Each lead is resolved independently through resolve_lead(); the rule enforced
    here is the only one a lead cannot know on its own, which is how many leads
    of its type are allowed at once.
    """
    per_lead = list(per_lead) or [{}]
    leads = [
        resolve_lead(cfg, entry, index, len(per_lead), overrides)
        for index, entry in enumerate(per_lead, start=1)
    ]

    cap = int(cfg.get("max_leads_per_type", DEFAULT_MAX_LEADS_PER_TYPE))
    counts = {}
    for lead in leads:
        counts[lead["type"]] = counts.get(lead["type"], 0) + 1
    over = sorted(t for t, c in counts.items() if c > cap)
    if over:
        raise ValueError(
            "that is %s, but max_leads_per_type is %d"
            % (" and ".join("%d %s leads" % (counts[t], t) for t in over), cap)
        )
    return leads


def resolve_lead(cfg, entry, index=1, of=1, overrides=None):
    """Flatten ONE lead into the numbers the geometry needs.

    Layering, later wins:

        lead_defaults.yaml `defaults`  ->  the type's own default block  ->
        this lead's own parameters (a panel tab, or the CLI flags)  ->
        `overrides`, applied to every lead of a set

    The type is settled first, across every layer, because it decides which
    default block sits at position two and which parameters are even legal.

    `level: T10` is turned into the z of that vertebral body here, so nothing
    downstream has to know what a level is; an explicit z_center in a LATER layer
    beats an earlier level, and vice versa, which is what lets the macro offer
    both controls without them fighting.

    Returns a plain dict, plus "index", "of", "label", and for an epidural lead
    "x" (absolute), for a DRG lead the ganglion "target".
    """
    entry = {k: v for k, v in (entry or {}).items() if v is not None}
    ovr = {k: v for k, v in (overrides or {}).items() if v is not None}

    kind = lead_type(cfg["defaults"], entry, ovr)
    type_defaults = cfg["drg_defaults"] if kind == "drg" else cfg["epidural_defaults"]
    layers = [cfg["defaults"], type_defaults, entry, ovr]

    params = {}
    for layer in layers:
        params.update(layer)
    params.pop("side", None)
    params["type"] = kind

    if kind == "drg":
        # Rejected if ASKED FOR, dropped if merely inherited. The distinction is
        # the point: `defaults:` describes what a lead is when nobody says
        # anything, and a DRG lead that quietly took its 8 contacts would be
        # lying about geometry it does not have. But a panel or a flag that asks
        # for six contacts on a DRG lead is asking for something this type
        # cannot do, and saying so is better than doing something else.
        asked = sorted(k for k in DRG_FIXED_PARAMS
                       if any(l.get(k) is not None for l in (type_defaults, entry, ovr)))
        if asked:
            raise ValueError(
                "a DRG lead is a copy of the lead RADO ships, not a swept one, "
                "so %s %s fixed by its meshes and cannot be set -- remove %s. "
                "(Dorsal and ventral leads still take all five.)"
                % (", ".join(asked), "is" if len(asked) == 1 else "are",
                   "it" if len(asked) == 1 else "them"))
        for k in DRG_FIXED_PARAMS:
            params.pop(k, None)

    allowed = COMMON_PARAMS | (DRG_PARAMS if kind == "drg" else EPIDURAL_PARAMS)
    unknown = [k for k in params if k not in PARAM_TYPES]
    if unknown:
        raise ValueError("unknown lead parameter(s) %s" % ", ".join(sorted(unknown)))
    wrong = [k for k in params if k not in allowed]
    if wrong:
        other = "a DRG lead" if kind != "drg" else "a dorsal or ventral lead"
        raise ValueError(
            "this is a %s lead, so %s mean%s nothing to it -- %s belong%s to %s"
            % (
                kind,
                ", ".join(sorted(wrong)),
                "" if len(wrong) > 1 else "s",
                "they" if len(wrong) > 1 else "it",
                "" if len(wrong) > 1 else "s",
                other,
            )
        )

    for k, v in list(params.items()):
        if isinstance(v, str) and PARAM_TYPES[k] is not str:
            params[k] = PARAM_TYPES[k](v)

    if kind == "drg":
        # Measured back off RADO's five STLs, once, and cached. Everything
        # downstream -- check_sane, describe, the panel's summary line, the
        # report -- reads these the same way it reads a dorsal lead's, so
        # nothing else has to learn that a DRG lead is different.
        params.update(rado_drg_hardware())

    if kind in EPIDURAL_TYPES:
        level = params.get("level")
        if level:
            if level not in cfg["levels"]:
                raise ValueError(
                    "no vertebral level named %r (have: %s)"
                    % (level, ", ".join(cfg["levels"]))
                )
            # A level named in the same or a later layer than an explicit
            # z_center wins; that is what "later wins" means when both are
            # present, and the only way to tell is which layer each came from.
            if _layer_of("z_center", layers) <= _layer_of("level", layers):
                params["z_center"] = float(cfg["levels"][level])
        params["x"] = MIDLINE_X + float(params["x_offset"])

    params["index"] = index
    params["of"] = of
    params.setdefault("label", default_label(params, index, of))
    check_sane(params)
    return params


def lead_type(*layers):
    """The lead's type, from the last layer that names one. Accepts `side:`.

    `side: dorsal | ventral` is what this repo said before DRG leads existed,
    and `--side` is still a CLI flag, so it is accepted as a spelling of `type`
    rather than becoming a confusing "unknown parameter" error. It cannot spell
    `drg`: a DRG lead is not a side of the canal.
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
                    "a side of the canal, write `type: drg`" % side
                )
            kind = side
    if kind not in TYPES:
        raise ValueError("type must be one of %s, not %r" % (", ".join(TYPES), kind))
    return kind


def default_label(params, index, of):
    """A short human name for one lead, used in the tree and on disk.

    "lead 1 (dorsal)", which the preview turns into "lead 1 (dorsal) — contact
    08" and "lead 1 (dorsal) — insulator". The same form whether there is one
    lead or four: the number is what ties a body in the tree to a tab in the
    panel and to a subdirectory under generated_leads/, and a lone lead that
    called itself something else would break that tie for no gain.
    """
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
            problems.append(
                "target must name a ganglion, like L3 or R2, not %r" % target
            )
    if problems:
        raise ValueError("; ".join(problems))


def array_length(params):
    """Length of the contact array, mm (first contact's start to last one's end)."""
    return (
        params["contacts"] * params["contact_length"]
        + (params["contacts"] - 1) * params["gap"]
    )


def total_length(params):
    return array_length(params) + 2 * params["tail"]


def describe(params):
    """One lead as a clinician would say it out loud."""
    kind = params["type"]
    common = "%d contacts, %.2f mm long, %.2f mm apart, %.2f mm diameter" % (
        params["contacts"],
        params["contact_length"],
        params["gap"],
        params["diameter"],
    )
    if kind == "drg":
        target = params["target"]
        side = "left" if target[0] == "L" else "right"
        shift = float(params.get("lateral_offset", 0.0))
        if abs(shift) < 1e-9:
            where = "centred on the ganglion"
        else:
            where = "%.1f mm %s the ganglion along the root" % (
                abs(shift),
                "distal to" if shift > 0 else "proximal to",
            )
        return ("RADO's own lead (%s), copied onto ganglion %s through the %s "
                "neuroforamen, %s") % (common, target, side, where)
    dx = float(params["x_offset"])
    if abs(dx) < 1e-9:
        where = "on the midline"
    else:
        where = "%.1f mm %s of midline" % (abs(dx), "left" if dx > 0 else "right")
    at = params.get("level") or "z = %.1f mm" % params["z_center"]
    return "%s, in the %s epidural space at %s, %s" % (common, kind, at, where)


def describe_set(leads):
    """The whole configuration in one line."""
    if len(leads) == 1:
        return describe(leads[0])
    counts = {}
    for lead in leads:
        counts[lead["type"]] = counts.get(lead["type"], 0) + 1
    return "%d leads: %s" % (
        len(leads),
        ", ".join("%d %s" % (counts[t], t) for t in TYPES if t in counts),
    )


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
    return {
        "epidural": Mesh.Mesh(os.path.join(stl_dir, EPIDURAL_STL)),
        "dura": Mesh.Mesh(os.path.join(stl_dir, DURA_STL)),
    }


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
            % os.path.relpath(path, REPO)
        )
    return json.load(open(path))


def corridor_for(params, foramen):
    """The measured corridor for this lead's target ganglion."""
    target = params["target"]
    if target not in foramen["ganglia"]:
        raise ValueError(
            "no ganglion named %r in this model (have: %s)"
            % (target, ", ".join(sorted(foramen["ganglia"])))
        )
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
            "not enough to follow it" % (len(stations), side, x, z_lo, z_hi)
        )
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
        epi, params["x"], side, z_lo, z_hi
    )
    diff = [
        abs(mc.polyval(stored, z) - mc.polyval(local_coef, z)) for z, _y, _t in stations
    ]
    info = {
        "stations": stations,
        "skipped": skipped,
        "max_disagreement_mm": max(diff) if diff else 0.0,
    }

    if abs(float(params["x_offset"])) < 1e-9:
        info["source"] = "stored midline centreline (epidural_corridor.json)"
        coef = stored
    else:
        info["source"] = "centreline re-measured at x = %.2f (%d slices)" % (
            params["x"],
            len(stations),
        )
        coef = local_coef
    info["coef"] = coef
    return (lambda z: mc.polyval(coef, z)), info


# -- foraminal (DRG) --------------------------------------------------------
#
# What survives here of the swept DRG lead is the MEASUREMENT, not the
# construction. drg_point_of(), drg_u_center(), ganglion_standoff() and
# wall_clearance() are gone with the sweep they served: a DRG lead is now a
# rigid copy of RADO's own geometry (rado_drg_lead.py), so there is no
# centreline to evaluate, no array centre to slide along a corridor, and no
# standoff to solve for -- the lead carries its own relationship to its
# ganglion.
#
# remeasure_foramen() stays, and earns it. It walks the foraminal block along
# the placed lead's own span and reports how much channel there is and whether
# the stored quartic still describes the block. Neither answer moves the lead
# any more; both are what explain one target being tighter than another, which
# is exactly the question eight copies of one lead raise.

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
            "the lead is off the end of the foramen" % (params["target"], u_lo, u_hi)
        )
    # Compared only over the range the quartic was FITTED to. Outside it the
    # polynomial is extrapolating into the block's taper, where the measured
    # centre swings by a millimetre or more over a millimetre of u -- a real
    # property of the taper, not a drift between the fit and the block, and
    # folding it into this number would make the check cry wolf on every lead
    # whose tip runs a little past the ganglion.
    lo, hi = corridor["usable_u"]
    fitted = [r for r in rows if lo - 1e-9 <= r[0] <= hi + 1e-9]
    diff = max(
        (
            max(abs(r[1] - mf.polyval(cy, r[0])), abs(r[2] - mf.polyval(cz, r[0])))
            for r in fitted
        ),
        default=0.0,
    )
    return rows, {
        "max_disagreement_mm": diff,
        "n_stations": len(rows),
        "n_extrapolated": len(rows) - len(fitted),
    }


def build_shapes(params, path_info):
    """The lead itself: [(kind, index, shape)], plus the segment plan.

    One line of geometry, and it is not written here -- make_scs_lead owns the
    sweep, so all three lead types differ only in the curve handed to it.
    path_info is (point_of, t_center) from whichever measurement applies.
    """
    _fc, _Mesh, _p, msl, _mc, _mf = _geometry_modules()
    point_of, t_center = path_info
    segments, t0, total = msl.segment_plan(
        params["contacts"],
        params["contact_length"],
        params["gap"],
        params["tail"],
        t_center,
    )
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


# Extra directions used to break ties in inside(). Deliberately not axis-aligned
# and not parallel to each other: an axis-aligned ray through a mesh built on a
# regular grid is far more likely to graze an edge than a skew one.
TIEBREAK_DIRECTIONS = (
    (0.0, -1.0, 0.0),
    (0.5773502692, 0.5773502692, 0.5773502692),
    (-0.5773502692, 0.5773502692, -0.5773502692),
)


def _parity_inside(mesh, point, direction):
    """One ray. Odd number of forward crossings => inside the material.

    foraminate() returns the intersections of the whole INFINITE line, not the
    forward ray, so the crossings behind the point have to be dropped by hand --
    counting all of them would make every point look outside.
    """
    hits = mesh.foraminate(point, direction)
    n = 0
    for v in hits.values():
        t = (
            (v[0] - point[0]) * direction[0]
            + (v[1] - point[1]) * direction[1]
            + (v[2] - point[2]) * direction[2]
        )
        if t > 1e-9:
            n += 1
    return n % 2 == 1


def inside(mesh, point, direction=(0.0, 1.0, 0.0)):
    """Is `point` inside the material of this mesh?

    RADO's compartments are hollow shells, so "inside the epidural body" means
    inside its fat, not inside its central cavity -- the parity test
    measure_corridor.py documents.

    WHY THIS VOTES INSTEAD OF TRUSTING ONE RAY. Crossing parity is exact in
    theory and brittle in practice: a ray that passes exactly through a shared
    triangle edge or vertex can register one crossing or two depending on
    floating-point luck, and that single miscount inverts the answer. It is rare
    but it does not stay rare when you cast thousands of rays -- an 8-contact
    dorsal lead samples 4284 surface points.

    It produced a real false failure. A lead with 0.97 mm of clearance, whose
    every other measurement was comfortable, failed containment on ONE point of
    4284 at (58.065, 86.321, 126.386). Seven of eight ray directions put that
    point inside the fat; only the default +Y disagreed, because it grazed an
    edge. The report then said "100.0% of points are inside; the other 0.0% are
    outside", which is what a rounded percentage does to a count of one.

    Escalation is ASYMMETRIC on purpose: a bare "inside" is returned
    immediately, and only "outside" is re-tested along TIEBREAK_DIRECTIONS.
    Two reasons. The failure mode observed is outside-when-actually-inside, a
    single spurious crossing; being wrongly called inside needs an even number
    of spurious crossings, which is much less likely. And the points being
    tested sit on a lead that is overwhelmingly in open fat, so nearly all of
    them take the fast path and the check stays affordable -- voting on every
    point would cost four times the ray casts for no benefit.
    """
    if _parity_inside(mesh, point, direction):
        return True
    votes = 1 + sum(1 for d in TIEBREAK_DIRECTIONS if _parity_inside(mesh, point, d))
    return votes * 2 > 1 + len(TIEBREAK_DIRECTIONS)


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
        if not (
            bb.XMin <= p[0] <= bb.XMax
            and bb.YMin <= p[1] <= bb.YMax
            and bb.ZMin <= p[2] <= bb.ZMax
        ):
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
    return best**0.5


# --------------------------------------------------------------------------
# validation -- one lead
# --------------------------------------------------------------------------
def validate_lead(
    params,
    meshes=None,
    corridor=None,
    foramen=None,
    max_samples=6000,
    stl_dir=DEFAULT_STL,
):
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
    numbers.update(
        lead_z=[z0, z1],
        contacts_z=[c0, c1],
        length_mm=length,
        array_mm=array_length(params),
        x=params["x"],
    )

    # -- 1. does the centreline reach this far? --------------------------
    usable = corridor["sides"][side]["usable_z"]
    numbers["usable_z"] = usable
    if c0 < usable[0] or c1 > usable[1]:
        checks.append(
            (
                "FAIL",
                "Reaches past the model",
                "The contacts would sit at z %.1f-%.1f mm, but this model's "
                "%s epidural space is only measured over z %.1f-%.1f mm. "
                "Part of the active lead would be outside the anatomy. Move "
                "the lead %.1f mm %s, or use fewer contacts."
                % (
                    c0,
                    c1,
                    side,
                    usable[0],
                    usable[1],
                    max(usable[0] - c0, c1 - usable[1]),
                    "rostral" if c0 < usable[0] else "caudal",
                ),
            )
        )
    elif z0 < usable[0] or z1 > usable[1]:
        checks.append(
            (
                "WARN",
                "Tail runs off the measured corridor",
                "Every contact is inside the measured anatomy (z %.1f-%.1f mm), "
                "but the passive tail reaches z %.1f-%.1f mm, past the "
                "measured range z %.1f-%.1f mm, so the last %.1f mm of "
                "insulator follows an extrapolated curve. Harmless for the "
                "field, worth knowing before meshing."
                % (
                    c0,
                    c1,
                    z0,
                    z1,
                    usable[0],
                    usable[1],
                    max(usable[0] - z0, z1 - usable[1]),
                ),
            )
        )
    else:
        checks.append(
            (
                "PASS",
                "Inside the modelled levels",
                "The lead spans z %.1f-%.1f mm, within the measured %s "
                "corridor z %.1f-%.1f mm." % (z0, z1, side, usable[0], usable[1]),
            )
        )

    # -- 2. build it, following the right curve ---------------------------
    try:
        y_of_z, cinfo = centreline_for(params, epi, corridor, (z0, z1))
    except ValueError as exc:
        checks.append(
            (
                "FAIL",
                "No usable channel here",
                "The %s epidural space could not be followed at this "
                "position: %s. This usually means the lead is off the end "
                "of the dural sac or beside a nerve-root sleeve rather "
                "than in the canal." % (side, exc),
            )
        )
        return {
            "ok": False,
            "checks": checks,
            "numbers": numbers,
            "shapes": None,
            "params": params,
        }

    numbers["centreline_source"] = cinfo["source"]
    numbers["centreline_vs_stored_max_mm"] = cinfo["max_disagreement_mm"]
    numbers["skipped_slices"] = len(cinfo["skipped"])
    shapes, segments, _z0, _total = build_shapes(
        params, (lambda z: (params["x"], y_of_z(z), z), params["z_center"])
    )

    if cinfo["skipped"]:
        checks.append(
            (
                "WARN",
                "Some slices could not be measured",
                "%d of %d sampled cross-sections along the lead did not give "
                "a clean channel (at z %s). Those are usually foramen levels "
                "where a nerve root leaves; the centreline is interpolated "
                "through them."
                % (
                    len(cinfo["skipped"]),
                    len(cinfo["skipped"]) + len(cinfo["stations"]),
                    ", ".join("%.0f" % z for z in cinfo["skipped"][:8]),
                ),
            )
        )

    # -- 3. does the lead fit the channel where it actually lies? ---------
    thicknesses = [t for _z, _y, t in cinfo["stations"]]
    t_min, t_med = min(thicknesses), sorted(thicknesses)[len(thicknesses) // 2]
    z_at_min = min(cinfo["stations"], key=lambda s: s[2])[0]
    clearance = t_min - params["diameter"]
    numbers.update(
        channel_min_mm=t_min,
        channel_median_mm=t_med,
        channel_min_at_z=z_at_min,
        clearance_mm=clearance,
        widest_lead_mm=t_min,
    )
    where = (
        "on the midline"
        if abs(float(params["x_offset"])) < 1e-9
        else "%.1f mm %s of midline"
        % (
            abs(float(params["x_offset"])),
            "left" if float(params["x_offset"]) > 0 else "right",
        )
    )
    if clearance < 0:
        checks.append(
            (
                "FAIL",
                "Too thick for this space",
                "The lead is %.2f mm across, but the %s epidural fat %s is "
                "only %.2f mm thick at its narrowest (z = %.0f mm). It would "
                "press through the dura. The widest lead that fits here is "
                "%.2f mm; the %s space is the roomier one at %.2f mm."
                % (
                    params["diameter"],
                    side,
                    where,
                    t_min,
                    z_at_min,
                    t_min,
                    "dorsal" if side == "ventral" else "ventral",
                    corridor["sides"]["dorsal" if side == "ventral" else "ventral"][
                        "thickness_min_mm"
                    ],
                ),
            )
        )
    elif clearance < TIGHT_CLEARANCE_MM:
        checks.append(
            (
                "WARN",
                "A tight fit",
                "The lead is %.2f mm across and the %s fat %s is %.2f mm at "
                "its narrowest (z = %.0f mm) -- only %.2f mm to spare, i.e. "
                "%.2f mm each side. It fits, but it is against both walls."
                % (
                    params["diameter"],
                    side,
                    where,
                    t_min,
                    z_at_min,
                    clearance,
                    0.5 * clearance,
                ),
            )
        )
    else:
        checks.append(
            (
                "PASS",
                "Fits the epidural space",
                "The %s fat %s is %.2f mm thick at its narrowest along this "
                "lead (median %.2f mm) against a %.2f mm lead: %.2f mm of "
                "clearance, %.2f mm each side."
                % (
                    side,
                    where,
                    t_min,
                    t_med,
                    params["diameter"],
                    clearance,
                    0.5 * clearance,
                ),
            )
        )

    # -- 4. is the built lead actually in the fat, and out of the dura? ---
    pts = sample_points(shapes, max_samples=max_samples)
    in_fat = [p for p in pts if inside(epi, p)]
    in_dura = [p for p in pts if inside(dura, p)]
    pct_fat = 100.0 * len(in_fat) / len(pts)
    pct_dura = 100.0 * len(in_dura) / len(pts)
    numbers.update(sampled_vertices=len(pts), pct_in_fat=pct_fat, pct_in_dura=pct_dura)

    # Report COUNTS, not percentages. "100.0% are inside; the other 0.0% are
    # outside" is what %.1f does to 1 stray point in 4284, and a verdict whose
    # own evidence reads as 100% pass is worse than no verdict at all.
    n_stray = len(pts) - len(in_fat)
    if pct_dura > 0:
        checks.append(
            (
                "FAIL",
                "Overlaps the dura",
                "%d of %d sampled points on the lead surface (%.2f%%) are "
                "inside the dura, between z %.1f and %.1f mm. The lead is not "
                "in the epidural space there -- it is in the thecal sac."
                % (
                    len(in_dura),
                    len(pts),
                    pct_dura,
                    min(p[2] for p in in_dura),
                    max(p[2] for p in in_dura),
                ),
            )
        )
    elif n_stray:
        stray = [p for p in pts if not inside(epi, p)]
        z_lo, z_hi = min(p[2] for p in stray), max(p[2] for p in stray)
        span = (
            "at z %.1f mm" % z_lo
            if z_hi - z_lo < 0.05
            else "between z %.1f and %.1f mm" % (z_lo, z_hi)
        )
        checks.append(
            (
                "FAIL",
                "Leaves the epidural fat",
                "%d of %d sampled points on the lead surface (%.2f%%) are "
                "outside the epidural fat, %s. The lead would be sitting "
                "partly in bone or in the foramen."
                % (n_stray, len(pts), 100.0 - pct_fat, span),
            )
        )
    else:
        checks.append(
            (
                "PASS",
                "Entirely within the epidural fat",
                "All %d sampled points on the lead surface are inside the "
                "epidural fat and none are in the dura." % len(pts),
            )
        )

    ok = not any(level == "FAIL" for level, _t, _m in checks)
    return {
        "ok": ok,
        "checks": checks,
        "numbers": numbers,
        "points": pts,
        "shapes": shapes,
        "segments": segments,
        "params": params,
    }


# Bone and disc. A lead tip that leaves the foraminal block has to be tested
# against these, because "outside the foramen" has two completely different
# meanings: out of the DISTAL MOUTH into paravertebral soft tissue, which is
# where a real DRG lead's tip goes, or INTO THE PEDICLE, which is not a
# placement. Named from ../ansys/tissue_map.yaml's `vertebra` and
# `intervertebral_disc` patterns.
BONE_STLS = ("T8-10 - V1-2.STL", "T8-10 - V1-3.STL", "T8-10 - V1-4.STL")
DISC_STLS = ("T8-10 - T_disk_9_10x-1.STL", "T8-10 - T_disk_10_11x-1.STL",
             "T8-10 - T_disk_11_12x-1.STL", "T8-10 - T_disk_12_Lx-1.STL")


def _sample_tris(tris, max_samples):
    """Deterministic vertex sample of one transformed body, deduplicated.

    An STL stores every facet's three vertices independently, so a body's raw
    vertex list repeats each corner five or six times. Deduplicating first
    means the sample covers the surface evenly instead of over-weighting
    whichever regions RADO tessellated finely, and it cuts the ray casts --
    the dominant cost of this whole check -- by about a factor of six.
    """
    import numpy as np
    v = np.unique(np.round(tris.reshape(-1, 3), 6), axis=0)
    if len(v) > max_samples:
        v = v[:: len(v) // max_samples + 1]
    return [tuple(float(c) for c in p) for p in v]


def _validate_drg(params, meshes, foramen, max_samples, stl_dir):
    """A DRG lead's checks. It is a COPY of RADO's lead, not a swept one.

    The questions are different from the swept version's, because two of them
    no longer need asking. Whether the contacts are the right size and spacing
    is not a question any more -- they are RADO's meshes and a rigid transform
    cannot change them. Neither is whether the trajectory matches RADO's: at the
    ganglion RADO built it against it IS RADO's, exactly, and everywhere else it
    is that same geometry carried over by a transform whose rigidity is checked
    arithmetically rather than measured geometrically.

    What is left is the part that genuinely differs from foramen to foramen:

        placement   the transform is rigid, its determinant is +1 or -1 as the
                    side demands, and every transformed body came out watertight
                    with OUTWARD normals. A cross-side copy is a reflection and
                    inverts normals unless the winding is reversed, so this is
                    checked rather than assumed, per body, by signed volume.
        clearance   nothing of the lead inside the ganglion it is stimulating,
                    and nothing inside the thecal sac. Either is a FAIL naming
                    which, because either is an injury rather than a placement.
        containment every sampled vertex inside the foraminal tissue -- and
                    where it is not, WHAT it is in. Out of the distal mouth into
                    unmodelled paravertebral space is what a real lead's tip
                    does and warns; into the pedicle or a disc fails.
        standoff    how far the lead surface sits off the ganglion surface,
                    measured point-to-triangle, against the same number for
                    RADO's own lead at its own ganglion. This is the number
                    "interfaces smoothly with the corresponding DRG" cashes out
                    as, and it is the reason the frame is built on the ganglion
                    rather than on the corridor: the eight cores are congruent
                    copies, so a same-side transform preserves it exactly.
        foramen     how much channel there is here, from measure_foramen.py's
                    stations, re-walked along the lead's own span. This no
                    longer DRIVES anything -- the geometry is copied, not swept
                    along the corridor -- but it is what explains a target being
                    tighter than another, so it is still measured and reported,
                    and the stored quartic is still checked against the block.
    """
    import numpy as np
    import rado_drg_lead as rdl

    checks, numbers = [], {}
    if foramen is None:
        foramen = load_foramen_corridors()
    try:
        corridor = corridor_for(params, foramen)
    except ValueError as exc:
        checks.append(("FAIL", "No such ganglion", str(exc)))
        return {"ok": False, "checks": checks, "numbers": numbers,
                "shapes": None, "tris": None, "params": params}

    target = params["target"]
    cache = meshes.setdefault("_rado_drg_cache", {})

    # -- 1. place it, and audit the transform ------------------------------
    placed = rdl.place(target, foramen, stl_dir,
                       float(params.get("lateral_offset", 0.0)),
                       float(params.get("y_offset", 0.0)),
                       float(params.get("z_offset", 0.0)), cache)
    hw = placed["hardware"]
    numbers.update(
        target=target,
        side=placed["side"],
        source_ganglion=placed["source"],
        mirrored=placed["mirrored"],
        transform_det=placed["det"],
        ganglion_z=corridor["ganglion_centroid"][2],
        ganglion_u=corridor["ganglion_u"],
        usable_u=corridor["usable_u"],
        hardware=hw,
        axis_vs_corridor_deg=placed["target_frame"]["axis_vs_corridor_deg"],
    )
    if placed["mirrored"]:
        numbers["sagittal"] = {
            "residual_rotation_deg": placed["sagittal"]["residual_rotation_deg"],
            "ganglion_asymmetry_mm": placed["sagittal"]["ganglion_asymmetry_mm"],
        }
    bad = [r for r in placed["orientation"] if not r["ok"]]
    if bad:
        checks.append((
            "FAIL", "The transformed geometry is not sound",
            "%d of the 5 copied bodies came out wrong: %s. A cross-side copy "
            "is a reflection, which inverts every facet normal unless the "
            "vertex winding is reversed with it, and an inside-out STL will "
            "not mesh. This is a bug in the transform, not something an offset "
            "can fix."
            % (len(bad), ", ".join(
                "%s (%s)" % (r["body"], "normals inward" if not r["normals_outward"]
                             else "not rigid" if not r["volume_preserved"]
                             else "%d open edges" % r["open_edges"]) for r in bad))))
        return {"ok": False, "checks": checks, "numbers": numbers,
                "shapes": None, "tris": None, "placed": placed, "params": params}
    checks.append((
        "PASS", "Copied onto the ganglion cleanly",
        "RADO's own 4-contact lead, carried from ganglion %s onto %s by a %s "
        "(determinant %+.0f). All five bodies came back watertight, the same "
        "volume and area to 1 part in a million, and with outward normals%s."
        % (placed["source"], target,
           "reflection through the sagittal plane" if placed["mirrored"]
           else "rotation within the same side", placed["det"],
           " after their windings were reversed" if placed["mirrored"] else "")))

    # -- 2. sample, by body, so a stray point can be attributed ------------
    per_body, pts = [], []
    share = max(400, max_samples // len(placed["bodies"]))
    for kind, idx, fn, tris in placed["bodies"]:
        p = _sample_tris(tris, share)
        per_body.append((kind, idx, fn, p))
        pts += p
    numbers["sampled_vertices"] = len(pts)

    block = mesh_named(meshes, corridor["block_stl"], stl_dir)
    drg = mesh_named(meshes, corridor["drg_core_stl"], stl_dir)
    dura = meshes["dura"]

    in_drg = inside_any(drg, pts)
    in_dura = inside_any(dura, pts)
    numbers["pct_in_ganglion"] = 100.0 * len(in_drg) / len(pts)
    numbers["pct_in_dura"] = 100.0 * len(in_dura) / len(pts)

    # -- 3. the two that are injuries rather than placements ---------------
    if in_drg:
        checks.append((
            "FAIL", "Passes through the ganglion",
            "%.2f%% of %d sampled points on the lead surface are inside "
            "ganglion %s itself. A DRG lead sits BESIDE its ganglion in the "
            "foraminal fat; one that goes through it is not a placement. "
            "Nudge it off with z_offset (rostral) or y_offset (dorsal)."
            % (numbers["pct_in_ganglion"], len(pts), target)))
    if in_dura:
        checks.append((
            "FAIL", "Enters the thecal sac",
            "%.2f%% of %d sampled points are inside the dura at %s. The lead "
            "has been pushed back into the canal rather than out through the "
            "foramen; reduce lateral_offset, or if it is zero this foramen "
            "cannot take the lead at its native depth."
            % (numbers["pct_in_dura"], len(pts), target)))
    if not in_drg and not in_dura:
        checks.append((
            "PASS", "Clear of the ganglion and the thecal sac",
            "None of the %d sampled points on the lead surface are inside "
            "ganglion %s, and none are inside the dura." % (len(pts), target)))

    # -- 4. containment, and what any stray point is actually in -----------
    in_block = set(inside_any(block, pts))
    numbers["pct_in_foraminal_tissue"] = 100.0 * len(in_block) / len(pts)
    by_kind = {}
    for kind, _idx, _fn, p in per_body:
        n_in = sum(1 for q in p if q in in_block)
        a, b = by_kind.get(kind, (0, 0))
        by_kind[kind] = (a + n_in, b + len(p))
    numbers["pct_in_tissue_by_kind"] = {
        k: 100.0 * a / b for k, (a, b) in by_kind.items()}
    stray = [q for q in pts if q not in in_block]

    if not stray:
        checks.append((
            "PASS", "Wholly inside the foraminal tissue",
            "All %d sampled points on the lead surface are inside the %s "
            "foraminal tissue." % (len(pts), target)))
    else:
        s = np.array(stray)
        u = np.abs(s[:, 0] - MIDLINE_X)
        depth = rdl.min_distance_to_surface(
            s, rdl.read_stl(os.path.join(stl_dir, corridor["block_stl"])))
        worst = max(rdl.min_distance_to_surface(s[i:i + 1], rdl.read_stl(
            os.path.join(stl_dir, corridor["block_stl"])))
            for i in range(0, len(s), max(1, len(s) // 20)))
        solid = []
        for fn in BONE_STLS + DISC_STLS:
            hit = inside_any(mesh_named(meshes, fn, stl_dir), stray)
            if hit:
                solid.append((fn, len(hit)))
        numbers.update(stray_points=len(stray), stray_u=[float(u.min()), float(u.max())],
                       stray_depth_mm=float(worst),
                       stray_in_bone=sum(n for _f, n in solid))
        if solid:
            checks.append((
                "FAIL", "Part of the lead is in bone",
                "%d of %d sampled points are outside the %s foraminal tissue "
                "AND inside %s. The lead is in the pedicle or a disc, not in "
                "the foramen."
                % (sum(n for _f, n in solid), len(pts), target,
                   ", ".join(f.replace("T8-10 - ", "").replace(".STL", "")
                             for f, _n in solid))))
        else:
            # The benign case, and the commonest one. Measured, not assumed:
            # every stray point is tested against all three vertebrae and all
            # four discs above before this branch is reached.
            checks.append((
                "WARN", "The tip protrudes from the mouth of the foramen",
                "%.1f%% of %d sampled points -- all of them at the distal end, "
                "%.1f to %.1f mm from the midline -- are outside the %s "
                "foraminal tissue, by at most %.2f mm. They are in none of the "
                "three vertebrae and none of the four discs: they are past the "
                "lateral mouth of the foramen, in space no RADO compartment "
                "occupies. A real DRG lead's tip does exactly this. It still "
                "matters for meshing, because that much of the lead will have "
                "no surrounding tissue to be meshed against."
                % (100.0 * len(stray) / len(pts), len(pts), u.min(), u.max(),
                   target, worst)))

    # -- 5. the standoff, which is what "interfaces smoothly" means --------
    core = rdl.read_stl(os.path.join(stl_dir, corridor["drg_core_stl"]))
    lead_v = rdl.vertices(placed)
    gap = rdl.min_distance_to_surface(lead_v[:: max(1, len(lead_v) // 8000)], core)
    native = rado_native_standoff(foramen, stl_dir, cache)
    numbers["ganglion_standoff_mm"] = gap
    numbers["rado_standoff_mm"] = native
    numbers["standoff_vs_rado_mm"] = gap - native
    checks.append((
        "PASS" if abs(gap - native) < 0.25 else "WARN",
        "Standoff from the ganglion",
        "The lead surface comes within %.3f mm of ganglion %s. RADO's own lead "
        "sits %.3f mm off ganglion %s, so this placement is %+.3f mm %s. %s"
        % (gap, target, native, placed["source"], gap - native,
           "further off" if gap > native else "closer",
           "The eight ganglion cores are congruent copies, so a same-side copy "
           "reproduces RADO's standoff exactly and a cross-side one differs "
           "only by the small difference between the left and right master "
           "bodies." if abs(gap - native) < 0.25 else
           "That is a bigger difference than the two master ganglion bodies "
           "account for; check the offsets.")))

    # -- 6. how much foramen there is here, and is the stored fit still true -
    u_lead = np.abs(lead_v[:, 0] - MIDLINE_X)
    span = (float(u_lead.min()), float(u_lead.max()))
    numbers["lead_u"] = list(span)
    try:
        rows, finfo = remeasure_foramen(block, corridor, params, span[0], span[1])
    except ValueError:
        rows, finfo = [], None
    if rows:
        tight = min(rows, key=lambda r: min(r[3], r[4]))
        numbers.update(
            foramen_min_section_mm=[float(tight[3]), float(tight[4])],
            foramen_min_at_u=float(tight[0]),
            centreline_vs_stored_max_mm=finfo["max_disagreement_mm"],
            extrapolated_stations=finfo["n_extrapolated"],
            centreline_source="stored quartic for %s (foraminal_corridors.json, "
                              "y rms %.3f mm, z rms %.3f mm)"
                              % (target, corridor["y_fit_rms_mm"],
                                 corridor["z_fit_rms_mm"]))
        room = min(tight[3], tight[4]) - hw["insulator_diameter"]
        numbers["channel_room_mm"] = float(room)
        if room < TIGHT_CLEARANCE_MM:
            checks.append((
                "WARN", "A tight foramen",
                "At %.1f mm from the midline the %s foramen is only %.2f mm by "
                "%.2f mm, against a %.2f mm sheath -- %.2f mm to spare. The "
                "lead is not being resized to fit; this is what the channel "
                "measures here."
                % (tight[0], target, tight[3], tight[4],
                   hw["insulator_diameter"], room)))

    # -- 7. root sleeves, which the foraminal block is not carved around ---
    lead_bb = (lead_v[:, 0].min(), lead_v[:, 0].max(), lead_v[:, 1].min(),
               lead_v[:, 1].max(), lead_v[:, 2].min(), lead_v[:, 2].max())
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
    if pierced:
        # Not a failure, and the reason is RADO's own lead: measured against
        # these same meshes at its own ganglion it clips three root bodies while
        # sitting wholly inside the foraminal tissue and wholly outside the
        # ganglion. RADO's compartments are not carved out of each other -- the
        # foraminal block CONTAINS the root and the ganglion -- so a small
        # overlap with a root sleeve is a property of the model rather than of
        # the placement. It still matters for meshing, so it is counted.
        checks.append((
            "WARN", "Grazes a nerve-root body",
            "%d sampled point%s of the lead surface fall inside %d nerve-root "
            "%s near this foramen (%s). RADO's own lead does the same at its "
            "own ganglion, because the foraminal block is not carved out around "
            "the roots -- but the mesher will have to resolve the overlap."
            % (sum(n for _f, n in pierced),
               "" if sum(n for _f, n in pierced) == 1 else "s", len(pierced),
               "body" if len(pierced) == 1 else "bodies",
               ", ".join(fn.replace("T8-10 - ", "")[:38] for fn, _n in pierced[:3]))))

    ok = not any(level == "FAIL" for level, _t, _m in checks)
    return {"ok": ok, "checks": checks, "numbers": numbers, "points": pts,
            "shapes": None, "tris": [(k, i, f, t) for k, i, f, t in placed["bodies"]],
            "placed": placed, "params": params}


def rado_native_standoff(foramen, stl_dir=DEFAULT_STL, cache=None):
    """How far RADO's own lead sits off its own ganglion. mm, measured once.

    The reference every other target's standoff is quoted against, and the one
    number in this module that says what "the RADO DRG lead was perfect" means
    quantitatively. Point-to-triangle, not vertex-to-vertex, because it is being
    compared against other targets at the 0.01 mm level.
    """
    import rado_drg_lead as rdl
    if cache is None:
        cache = {}
    if "native_standoff" not in cache:
        info = foramen["ganglia"][rdl.NATIVE_TARGET]
        core = rdl.read_stl(os.path.join(stl_dir, info["drg_core_stl"]))
        lead = rdl.vertices({"bodies": rdl.load_rado_bodies(stl_dir, cache)})
        cache["native_standoff"] = rdl.min_distance_to_surface(
            lead[:: max(1, len(lead) // 8000)], core)
    return cache["native_standoff"]


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


def lead_mesh(result, stl_dir=DEFAULT_STL):
    """One copied lead's five bodies as a single FreeCAD Mesh. For collisions.

    A DRG lead has no Part shape to fuse -- it is transformed triangles -- so
    the pairwise check gets a mesh instead. Mesh.addMesh() concatenates rather
    than booleans, which is all that is wanted here: the five bodies do not
    overlap each other and the parity test does not care whether they are one
    surface or five.
    """
    _fc, Mesh, _p, _msl, _mc, _mf = _geometry_modules()
    out = Mesh.Mesh()
    for _kind, _idx, _fn, tris in result["tris"]:
        out.addMesh(Mesh.Mesh([tuple(tuple(float(c) for c in v) for v in t)
                               for t in tris]))
    return out


def pair_clearance(a_result, b_result):
    """How close two built leads come, and whether they overlap. A dict.

    TWO MEASUREMENTS, and which one is used depends on what the leads are.

    Between two SWEPT leads the gap is the exact shape-to-shape distance, not a
    sampled one, because "do these two 1.3 mm cylinders touch" is a question OCC
    answers exactly and a vertex sample answers only probably. When the gap is
    zero the overlap is quantified two ways: the shared volume, and the fraction
    of one lead's sampled surface inside the other, which is the number that
    reads clinically ("a third of lead 2 is inside lead 1").

    When either lead is a COPIED DRG lead there is no solid to hand OCC -- it is
    35000 triangles, and makeShapeFromMesh on that is minutes of work to answer
    a question about two objects that are usually centimetres apart. So the gap
    is measured point-to-triangle between the two surfaces instead, which is
    exact for the nearest approach up to the tessellation the meshes already
    are, and the overlap is measured by the same crossing-parity test the
    containment checks use. `gap_exact` records which was done, so a report
    never quotes a sampled number as though OCC had said it.

    The only thing lost is `overlap_mm3` for a mesh pair: two intersecting
    meshes have no cheap shared volume here, so it comes back NaN and the
    percentage carries the verdict. A pair of leads that intersect at all is a
    FAIL either way; how many cubic millimetres they share changes nothing.
    """
    import numpy as np
    import rado_drg_lead as rdl

    if is_copy(a_result) or is_copy(b_result):
        out = {"gap_mm": 0.0, "overlap_mm3": float("nan"), "pct_a_in_b": 0.0,
               "gap_exact": False}
        a_pts = (np.array(a_result["points"]) if a_result.get("points")
                 else np.array(sample_points(a_result["shapes"])))
        if is_copy(b_result):
            b_tris = np.vstack([t for _k, _i, _f, t in b_result["tris"]])
            out["gap_mm"] = rdl.min_distance_to_surface(a_pts, b_tris)
            mesh_b = lead_mesh(b_result)
        else:
            b_pts = (np.array(b_result["points"]) if b_result.get("points")
                     else np.array(sample_points(b_result["shapes"])))
            a_tris = np.vstack([t for _k, _i, _f, t in a_result["tris"]])
            out["gap_mm"] = rdl.min_distance_to_surface(b_pts, a_tris)
            mesh_b = None
        if out["gap_mm"] > 1e-7:
            return out
        if mesh_b is None:
            mesh_b = lead_mesh(b_result) if is_copy(b_result) else None
        if mesh_b is not None:
            n_in = len(inside_any(mesh_b, [tuple(p) for p in a_pts]))
            out["pct_a_in_b"] = 100.0 * n_in / len(a_pts) if len(a_pts) else 0.0
        return out

    solid_a = lead_solid(a_result["shapes"])
    solid_b = lead_solid(b_result["shapes"])
    gap = solid_a.distToShape(solid_b)[0]
    out = {"gap_mm": gap, "overlap_mm3": 0.0, "pct_a_in_b": 0.0, "gap_exact": True}
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


def validate_set(
    leads,
    meshes=None,
    corridor=None,
    foramen=None,
    max_samples=6000,
    stl_dir=DEFAULT_STL,
):
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

    results = [
        validate_lead(
            lead,
            meshes=meshes,
            corridor=corridor,
            foramen=foramen,
            max_samples=max_samples,
            stl_dir=stl_dir,
        )
        for lead in leads
    ]

    checks, pairs = [], []
    built = [(i, r) for i, r in enumerate(results) if is_built(r)]
    if len(built) < 2:
        if len(leads) > 1:
            checks.append(
                (
                    "WARN",
                    "Leads could not be compared",
                    "Fewer than two of the %d leads could be built, so they "
                    "have not been checked against each other." % len(leads),
                )
            )
    else:
        for ai in range(len(built)):
            for bi in range(ai + 1, len(built)):
                i, a = built[ai]
                j, b = built[bi]
                info = pair_clearance(a, b)
                info.update(
                    a=i, b=j, a_label=leads[i]["label"], b_label=leads[j]["label"]
                )
                pairs.append(info)

        worst = min(pairs, key=lambda p: p["gap_mm"])
        touching = [p for p in pairs if p["gap_mm"] <= 1e-7]
        tight = [p for p in pairs if 1e-7 < p["gap_mm"] < TIGHT_LEAD_GAP_MM]
        if touching:
            first = touching[0]
            checks.append(
                (
                    "FAIL",
                    "Two leads occupy the same space",
                    "%s and %s intersect: %s%.1f%% of the first lead's sampled "
                    "surface is inside the second. Two leads cannot both be there. "
                    "Move them apart -- for a dorsal pair, increase the difference "
                    "between their x_offset values; for a DRG lead, change its "
                    "target ganglion or its offsets."
                    % (
                        first["a_label"],
                        first["b_label"],
                        ("they share %.3f mm3 of volume, and "
                         % first["overlap_mm3"])
                        if first["overlap_mm3"] == first["overlap_mm3"] else "",
                        first["pct_a_in_b"],
                    ),
                )
            )
        elif tight:
            first = min(tight, key=lambda p: p["gap_mm"])
            checks.append(
                (
                    "WARN",
                    "Two leads are very close",
                    "%s and %s come within %.2f mm of each other. They do not touch, but "
                    "a mesher given a gap that small between two curved conductors will "
                    "produce sliver elements there, and the field between them will be "
                    "dominated by the gap rather than by the anatomy."
                    % (first["a_label"], first["b_label"], first["gap_mm"]),
                )
            )
        else:
            checks.append(
                (
                    "PASS",
                    "The leads clear each other",
                    "The closest approach between any two of the %d leads is %.2f mm "
                    "(%s and %s). None of them intersect."
                    % (len(leads), worst["gap_mm"], worst["a_label"], worst["b_label"]),
                )
            )

    ok = all(r["ok"] for r in results) and not any(
        level == "FAIL" for level, _t, _m in checks
    )
    return {"ok": ok, "leads": results, "pairs": pairs, "checks": checks}


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------
def report_lines(params, result):
    """One lead's validation as text, for the report file and the macro's panel."""
    n = result["numbers"]
    lines = [
        "%s -- %s" % (params["label"], "PASS" if result["ok"] else "FAIL"),
        describe(params),
        "",
    ]
    if params.get("note"):
        lines += [params["note"], ""]
    for level, title, message in result["checks"]:
        lines.append("[%-4s] %s" % (level, title))
        for chunk in wrap(message, 74):
            lines.append("        " + chunk)
    lines.append("")
    lines.append("measurements")
    if params["type"] == "drg":
        hw = n.get("hardware") or {}
        lines.append(
            "   source      RADO's own lead at ganglion %s, carried here by a %s"
            % (n.get("source_ganglion", "?"),
               "REFLECTION through the sagittal plane (det %+.0f)"
               % n.get("transform_det", 0) if n.get("mirrored")
               else "rotation within the same side (det %+.0f)"
               % n.get("transform_det", 0)))
        if n.get("sagittal"):
            lines.append(
                "               residual after the mirror: %.2f deg of rotation; "
                "the counterpart ganglia differ by %.2f mm"
                % (n["sagittal"]["residual_rotation_deg"],
                   n["sagittal"]["ganglion_asymmetry_mm"]))
        if hw:
            lines.append(
                "   hardware    %d contacts, %.3f mm long, %.3f mm pitch, %.3f mm "
                "across; sheath %.3f mm; FIXED, measured off RADO's meshes"
                % (hw["contacts"], hw["contact_length"], hw["pitch"],
                   hw["diameter"], hw.get("insulator_diameter", float("nan"))))
            lines.append(
                "               electrode numbers run %s medially to laterally "
                "-- RADO's numbering is not spatial"
                % " -> ".join(str(i) for i in hw["spatial_order"]))
        if "lead_u" in n:
            lines.append(
                "   lead        %.1f .. %.1f mm from the midline"
                % (n["lead_u"][0], n["lead_u"][1]))
            lines.append(
                "   target      ganglion %s (%s side), spanning %.1f .. %.1f mm "
                "from the midline at z = %.1f"
                % (n["target"], n["side"], n["ganglion_u"][0], n["ganglion_u"][1],
                   n["ganglion_z"]))
        if "pct_in_foraminal_tissue" in n:
            lines.append(
                "   containment %.2f%% of %d sampled vertices in foraminal tissue, "
                "%.2f%% in the ganglion, %.2f%% in dura"
                % (n["pct_in_foraminal_tissue"], n["sampled_vertices"],
                   n["pct_in_ganglion"], n["pct_in_dura"]))
            by = n.get("pct_in_tissue_by_kind") or {}
            if by:
                lines.append(
                    "               by body: %s"
                    % ", ".join("%s %.2f%%" % (k, v) for k, v in sorted(by.items())))
        if "stray_points" in n:
            lines.append(
                "   outside     %d vertices, %.1f .. %.1f mm from the midline, at "
                "most %.2f mm beyond the foraminal tissue; %d of them in bone or disc"
                % (n["stray_points"], n["stray_u"][0], n["stray_u"][1],
                   n["stray_depth_mm"], n["stray_in_bone"]))
        if "ganglion_standoff_mm" in n:
            lines.append(
                "   standoff    %.3f mm from the lead surface to ganglion %s; "
                "RADO's own is %.3f mm (%+.3f mm)"
                % (n["ganglion_standoff_mm"], n["target"], n["rado_standoff_mm"],
                   n["standoff_vs_rado_mm"]))
        if "foramen_min_section_mm" in n:
            lines.append(
                "   foramen     narrowest %.2f x %.2f mm at %.1f mm out, against a "
                "%.3f mm sheath: %.2f mm to spare"
                % (n["foramen_min_section_mm"][0], n["foramen_min_section_mm"][1],
                   n["foramen_min_at_u"], hw.get("insulator_diameter", float("nan")),
                   n.get("channel_room_mm", float("nan"))))
        if "axis_vs_corridor_deg" in n:
            lines.append(
                "   frame       the ganglion's long axis and the fitted foraminal "
                "corridor's tangent agree to %.2f deg here"
                % n["axis_vs_corridor_deg"])
    else:
        lines.append(
            "   lead        z %.1f .. %.1f mm (%.1f mm), contacts z %.1f .. %.1f"
            % (
                n["lead_z"][0],
                n["lead_z"][1],
                n["length_mm"],
                n["contacts_z"][0],
                n["contacts_z"][1],
            )
        )
        lines.append(
            "   position    x = %.2f mm (midline %.2f%+.2f)"
            % (n["x"], MIDLINE_X, float(params["x_offset"]))
        )
        if "channel_min_mm" in n:
            lines.append(
                "   channel     min %.2f mm at z %.1f, median %.2f mm; lead %.2f mm; "
                "clearance %+.2f mm"
                % (
                    n["channel_min_mm"],
                    n["channel_min_at_z"],
                    n["channel_median_mm"],
                    params["diameter"],
                    n["clearance_mm"],
                )
            )
            lines.append(
                "   containment %.1f%% of %d sampled vertices in epidural fat, "
                "%.1f%% in dura"
                % (n["pct_in_fat"], n["sampled_vertices"], n["pct_in_dura"])
            )
    if "centreline_source" in n:
        lines.append("   centreline  %s" % n["centreline_source"])
        lines.append(
            "   local measurement vs stored fit: max %.3f mm apart"
            % n["centreline_vs_stored_max_mm"]
        )
    return lines


def report_set_lines(name, leads, outcome):
    """The whole configuration's validation as text."""
    lines = [
        "%s -- %s" % (name, "PASS" if outcome["ok"] else "FAIL"),
        describe_set(leads),
        "",
    ]
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
            lines.append(
                "   %-22s %-22s %7.2f mm%s"
                % (
                    p["a_label"],
                    p["b_label"],
                    p["gap_mm"],
                    "   INTERSECTING" if p["gap_mm"] <= 1e-7 else "",
                )
            )
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
    """Measure a placed DRG lead against RADO's own, body by body. Lines.

    WHAT THIS CHECK IS NOW FOR. It used to be the DRG trajectory's only
    independent validation: a lead swept along the ray-cast foraminal corridor
    was compared against the one RADO placed by hand, and agreement meant the
    corridor was right. There is no swept trajectory left to validate -- the
    lead IS RADO's -- so the question changes to the one a copy raises instead:

        did the transform CHANGE the geometry, or only move it?

    A rigid map must preserve every internal measurement exactly. So this
    compares the placed lead against RADO's originals on the quantities a rigid
    map cannot change -- each body's signed volume, its surface area, and the
    full set of contact-to-contact centre distances -- and on the one it must
    change, position. At the native ganglion with no offsets the transform is
    the identity and EVERY displacement must be exactly zero; that is the
    strongest single assertion available about this whole mechanism, and it is
    what test_regression.py's DRG case turns into a pass or a fail.

    At any other ganglion the displacements are simply where the lead went, and
    the internal measurements carry the check on their own.
    """
    import numpy as np
    import rado_drg_lead as rdl

    lines = []
    if params["type"] != "drg":
        return ["--compare-rado only means anything for a DRG lead."]
    if not is_copy(result):
        return ["this DRG lead has no copied geometry to compare."]

    placed = result["placed"]
    native = rdl.load_rado_bodies(stl_dir)
    by_file = {fn: tris for _k, _i, fn, tris in native}

    lines += [
        "the COPY against RADO's originals in STL_files/",
        "(a rigid transform cannot change a volume, an area or an internal "
        "distance; it can only change a position)",
        "",
        "   %-26s %10s %10s %10s %10s %14s"
        % ("body", "vol orig", "vol copy", "area orig", "area copy",
           "centroid moved"),
    ]
    worst_v = worst_a = 0.0
    moves = []
    for _kind, _idx, fn, tris in result["tris"]:
        src = by_file[fn]
        v0, v1 = rdl.signed_volume(src), rdl.signed_volume(tris)
        a0, a1 = rdl.surface_area(src), rdl.surface_area(tris)
        move = float(np.linalg.norm(rdl.area_centroid(tris) - rdl.area_centroid(src)))
        moves.append(move)
        worst_v = max(worst_v, abs(v1 - v0))
        worst_a = max(worst_a, abs(a1 - a0))
        lines.append("   %-26s %10.5f %10.5f %10.4f %10.4f %11.4f mm"
                     % (fn, v0, v1, a0, a1, move))
    lines += [
        "",
        "   worst volume change %.2e mm3, worst area change %.2e mm2 -- both "
        "must be zero to rounding for the map to have been rigid"
        % (worst_v, worst_a),
    ]

    # The internal geometry: every contact-to-contact centre distance. Six
    # numbers for four contacts, and a rigid map preserves all six. Note what
    # this does NOT catch on its own -- a reflection preserves them too, and
    # would leave a mirror-image lead looking identical here. That is precisely
    # why the placement check audits the winding separately: the two together
    # pin the geometry, neither alone does.
    def spread(bodies):
        cs = [rdl.area_centroid(t) for k, _i, _f, t in bodies if k == "contact"]
        return sorted(float(np.linalg.norm(cs[i] - cs[j]))
                      for i in range(len(cs)) for j in range(i + 1, len(cs)))

    d0, d1 = spread(native), spread(result["tris"])
    lines += [
        "   contact-to-contact centre distances, sorted:",
        "      RADO   %s" % "  ".join("%.4f" % d for d in d0),
        "      copy   %s" % "  ".join("%.4f" % d for d in d1),
        "      worst difference %.2e mm" % max(abs(a - b) for a, b in zip(d0, d1)),
        "",
    ]

    if placed["target"] == placed["source"] and not any(
            abs(v) > 1e-12 for v in placed["offsets"].values()):
        ok = max(moves) == 0.0 and worst_v == 0.0 and worst_a == 0.0
        lines.append(
            "   this is the NATIVE placement (%s, no offsets), so the transform "
            "is the identity and every body must be exactly where RADO left it: "
            "%s -- largest centroid displacement %.2e mm"
            % (placed["source"], "CONFIRMED" if ok else "FAILED", max(moves)))
    else:
        lines.append(
            "   the lead has been carried from %s to %s%s, so the displacements "
            "above are the placement; the volumes, areas and internal distances "
            "are what had to stay put."
            % (placed["source"], placed["target"],
               " through the sagittal mirror" if placed["mirrored"] else ""))
    return lines


# --------------------------------------------------------------------------
# outputs: STLs on disk, and a disposable preview in a live document
# --------------------------------------------------------------------------
def lead_dir_name(params):
    """The subdirectory one lead of a multi-lead set is exported into."""
    if params["type"] == "drg":
        return "lead%d_drg_%s" % (params["index"], params["target"])
    return "lead%d_%s" % (params["index"], params["type"])


def default_tag(params):
    """A filing name for a single lead, when nobody supplied one.

    Only used by callers that want one; --tag and the panel's export dialog
    still win. For a DRG lead it names the ganglion, because that is now the
    ONLY thing that distinguishes one DRG lead from another -- the hardware is
    RADO's and identical on every target.
    """
    if params["type"] == "drg":
        return "drg_%s_rado" % params["target"]
    return "%s_z%.0f_%dc" % (params["type"], params["z_center"], params["contacts"])


def export_set(name, leads, results, out_root=DEFAULT_OUT):
    """Write the STLs for a whole set of leads. Returns (dir, [relative files]).

    `name` is just the directory: whatever the user typed in the panel's export
    dialog or passed as --tag. It is a filing decision, not a configuration --
    nothing reads it back.

    Layout, and it is deliberately not uniform:

        one lead   generated_leads/<name>/SCS Lead Electrode 1.stl, ...
        several    generated_leads/<name>/lead1_dorsal/SCS Lead Electrode 1.stl, ...

    A single lead keeps the flat layout it has always had, because
    generated_leads/dorsal_z110_8c/ and ventral_z110_8c/ are committed at those
    exact paths and are what test_regression.py compares against byte for byte;
    moving them would make that test untestable. A set gets one subdirectory per
    lead, plus a leads.txt manifest, because the FILENAMES cannot carry the
    distinction: they are RADO's own -- "SCS Lead Electrode 1.stl", "SCS Lead
    Insulator.stl" -- and ../ansys/tissue_map.yaml matches the object Names
    those sanitise to, so renaming them costs the contacts their silver and the
    insulator its material. The directory carries it instead.

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
    root = os.path.join(out_root, name)
    if not os.path.isdir(root):
        os.makedirs(root)
    several = len(leads) > 1

    doc = FreeCAD.newDocument("scs_lead_export")
    written = []
    manifest = ["%s -- %s" % (name, describe_set(leads)), ""]
    try:
        for params, result in zip(leads, results):
            outdir = os.path.join(root, lead_dir_name(params)) if several else root
            if not os.path.isdir(outdir):
                os.makedirs(outdir)
            manifest += [
                "%s  ->  %s" % (params["label"], os.path.relpath(outdir, root)),
                "   %s" % describe(params),
            ]
            if is_copy(result):
                # A copied DRG lead is written by rado_drg_lead's own binary
                # writer, not through FreeCAD. Not an optimisation: Mesh.export
                # would re-tessellate and re-derive normals, and the whole
                # claim about this geometry is that it is RADO's facets moved
                # rigidly and nothing else. Writing the arrays out preserves
                # that, and the winding fix travels with them.
                import rado_drg_lead as rdl
                for kind, idx, fn, tris in result["tris"]:
                    path = os.path.join(outdir, fn)
                    rdl.write_stl(path, tris)
                    written.append(os.path.relpath(path, root))
                hw = result["numbers"]["hardware"]
                manifest += [
                    "   copied from RADO's own lead at ganglion %s%s"
                    % (result["placed"]["source"],
                       ", mirrored through the sagittal plane"
                       if result["placed"]["mirrored"] else ""),
                    "   RADO's electrode numbers run %s medially to laterally --"
                    % " -> ".join(str(i) for i in hw["spatial_order"]),
                    "   they are NOT in file order, so pick a bipolar pair by "
                    "that sequence, not by the filenames.",
                ]
                manifest.append("")
                continue
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
    # A manifest for a SET, as always -- and also for any set containing a
    # copied DRG lead, however few leads that is, because a copied lead's
    # filenames carry RADO's own electrode numbering and that numbering is not
    # spatial. Someone picking contacts 1 and 2 for a bipolar pair off the
    # filenames would get the second and third along the root, not the first
    # two. The manifest is where that is written down.
    if several or any(is_copy(r) for r in results):
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


def preview_set(doc, leads, results, colours=None):
    """Put a whole set of leads into a live document as a disposable preview.

    They all land in ONE group, labelled "13 SCS Leads" -- every lead's contacts
    and its insulator, however many leads there are. One group rather than one
    per lead because that is what the tree is for: a lead is hardware, all the
    hardware is in one place, and the lead number lives in each body's Label
    ("lead 2 (drg) — contact 03") where it can be read without expanding
    anything. The Label does not say "preview"; what is in the tree is leads.

    Replaces any previous preview, so trying ten arrangements does not silt the
    tree up. The bodies are Part::Feature -- parametric shapes, a few kB -- not
    imported meshes, so a preview costs nothing and leaves nothing on disk. Which
    is the whole point: the STLs are the artefact, the document is not.

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
    FreeCAD, Mesh, _p, msl, _mc, _mf = _geometry_modules()
    if not FreeCAD.GuiUp and getattr(doc, "FileName", ""):
        raise RuntimeError(
            "refusing to preview into %s with no GUI running: a headless save "
            "would strip the document's colours. Preview from a running FreeCAD, "
            "or export STLs instead." % os.path.basename(doc.FileName)
        )

    clear_preview(doc)
    group = doc.addObject("App::DocumentObjectGroup", PREVIEW_GROUP)
    group.Label = PREVIEW_GROUP_LABEL

    made = []
    for params, result in zip(leads, results):
        if not is_built(result):
            continue
        n = params["index"]
        if is_copy(result):
            # A copied DRG lead previews as Mesh::Feature, because that is what
            # it is. Converting the 35000 triangles to a Part shape to make the
            # tree uniform would cost minutes and change the geometry into
            # something that is no longer RADO's facets. The contact NUMBER kept
            # here is RADO's own file number, not a spatial rank -- see the
            # export manifest -- so the label says so rather than implying an
            # order the numbers do not carry.
            for kind, idx, _fn, tris in result["tris"]:
                obj = doc.addObject(
                    "Mesh::Feature",
                    "SCS_Preview_Lead%d_%s_%02d"
                    % (n, "Contact" if kind == "contact" else "Insulator", idx))
                obj.Mesh = Mesh.Mesh([tuple(tuple(float(c) for c in v) for v in t)
                                      for t in tris])
                obj.Label = ("%s — contact %02d (RADO's numbering)" % (params["label"], idx)
                             if kind == "contact"
                             else "%s — insulator" % params["label"])
                made.append(("electrode_contact" if kind == "contact"
                             else "lead_insulation", obj))
            continue
        for kind, idx, shape in result["shapes"]:
            if kind != "contact":
                continue
            obj = doc.addObject(
                "Part::Feature", "SCS_Preview_Lead%d_Contact_%02d" % (n, idx)
            )
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
    point. Cleared with clear_rado_lead(). There is no button for it in the
    panel; it is two lines in FreeCAD's Python console, and drop_rado_lead.py's
    docstring has them.
    """
    FreeCAD, Mesh, _p, _msl, _mc, _mf = _geometry_modules()
    if not FreeCAD.GuiUp and getattr(doc, "FileName", ""):
        raise RuntimeError(
            "refusing to import into %s with no GUI running: a headless save "
            "would strip the document's colours." % os.path.basename(doc.FileName)
        )
    clear_rado_lead(doc)
    before = len(doc.Objects)
    for fn in RADO_LEAD_STLS:
        Mesh.insert(os.path.join(stl_dir, fn), doc.Name)
    added = doc.Objects[before:]
    group = doc.addObject("App::DocumentObjectGroup", RADO_REFERENCE_GROUP)
    group.Label = "REFERENCE: RADO's own 4-contact DRG lead"
    for fn, obj in zip(RADO_LEAD_STLS, added):
        obj.Label = (
            "RADO DRG lead insulator"
            if "Insulator" in fn
            else "RADO DRG lead contact %s" % os.path.splitext(fn)[0][-1]
        )
        if colours and obj.ViewObject is not None:
            rgb = colours.get(
                "lead_insulation" if "Insulator" in fn else "electrode_contact"
            )
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

    The lead tissues stay in that map even though the anatomy document holds no
    lead: a previewed lead needs its silver and its black, and an exported STL
    needs its conductivity when it reaches Ansys. It is the DOCUMENT that has no
    lead in it, not the model.

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
        return {
            t["name"]: t["color_rgb"]
            for t in tissues
            if t["name"] in ("electrode_contact", "lead_insulation") and t["color_rgb"]
        }
    except Exception:
        return {}


# --------------------------------------------------------------------------
def main():
    """One lead, from flags. No names, no catalogue -- see the module docstring."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument(
        "--defaults",
        default=DEFAULTS_YAML,
        help="starting values; default lead_defaults.yaml",
    )
    ap.add_argument("--corridor", default=DEFAULT_CORRIDOR)
    ap.add_argument("--foramen", default=DEFAULT_FORAMEN)
    ap.add_argument("--stl-dir", default=DEFAULT_STL)
    ap.add_argument("--out", default=DEFAULT_OUT)
    ap.add_argument("--report", default=DEFAULT_REPORT)
    ap.add_argument(
        "--show",
        action="store_true",
        help="resolve to numbers and stop; no FreeCAD needed",
    )
    ap.add_argument("--validate", action="store_true", help="measure the fit")
    ap.add_argument(
        "--compare-rado",
        action="store_true",
        help="measure a DRG lead against RADO's own 4-contact lead",
    )
    ap.add_argument(
        "--export", action="store_true", help="validate, then write STLs if it passes"
    )
    ap.add_argument(
        "--force",
        action="store_true",
        help="export even if validation failed (it will not fit)",
    )
    ap.add_argument(
        "--tag",
        default="custom_lead",
        help="output directory name under generated_leads/",
    )
    ap.add_argument("--max-samples", type=int, default=6000)
    for flag, kind in (
        ("--type", str),
        ("--side", str),
        ("--contacts", int),
        ("--contact-length", float),
        ("--gap", float),
        ("--diameter", float),
        ("--x-offset", float),
        ("--z-center", float),
        ("--tail", float),
        ("--level", str),
        ("--target", str),
        ("--lateral-offset", float),
        ("--y-offset", float),
        ("--z-offset", float),
    ):
        ap.add_argument(
            flag,
            type=kind,
            default=None,
            help="lead parameter; unset means lead_defaults.yaml",
        )
    args = ap.parse_args(script_args())

    out = []

    def say(fmt, *a):
        out.append(fmt % a if a else fmt)

    cfg = parse_defaults_file(args.defaults)
    params = {
        k: getattr(args, k.replace("-", "_"))
        for k in (
            "type",
            "side",
            "contacts",
            "contact_length",
            "gap",
            "diameter",
            "x_offset",
            "z_center",
            "tail",
            "level",
            "target",
            "lateral_offset",
            "y_offset",
            "z_offset",
        )
    }

    if not (args.show or args.validate or args.export or args.compare_rado):
        say("Nothing to do. Say what you want done with the lead:")
        say("")
        say("    --show           resolve the parameters to numbers and stop")
        say("    --validate       measure whether it fits the anatomy")
        say("    --export         validate, then write its STLs")
        say("    --compare-rado   measure a DRG lead against RADO's own")
        say("")
        say("A DORSAL or VENTRAL lead is flags: --type, --contacts,")
        say("--contact-length, --gap, --diameter, --tail, and then")
        say("--level/--z-center/--x-offset. A DRG lead takes NONE of the")
        say("hardware flags -- it is a copy of the lead RADO ships and its")
        say("contacts are fixed -- only --target and the three nudges")
        say("--lateral-offset/--y-offset/--z-offset.")
        say("")
        say("Unchanged:")
        say(
            "one. Anything not given comes from %s.",
            os.path.relpath(args.defaults, REPO),
        )
        write_report(args.report, out)
        return 2

    # Fail here rather than three functions deep with an ImportError traceback:
    # measuring anything needs FreeCAD, and the fix is a different interpreter,
    # not a different lead.
    if args.validate or args.export or args.compare_rado:
        try:
            _geometry_modules()
        except ImportError as exc:
            say(
                "Measuring a lead needs FreeCAD, and this interpreter has none (%s).",
                exc,
            )
            say("")
            say("Run it under freecadcmd instead:")
            say(
                '    BUILD_LEAD_CONFIG_ARGS="%s" freecadcmd src/freecad/%s',
                " ".join(script_args()),
                os.path.basename(__file__),
            )
            say("")
            say("--show needs no FreeCAD and works here.")
            say("")
            say("RESULT: aborted, no FreeCAD in this interpreter")
            write_report(args.report, out)
            return 2

    try:
        leads = resolve_leads(cfg, [params])
    except (KeyError, ValueError) as exc:
        # An impossible lead is a user mistake, not a bug: say what is wrong.
        say("=" * 72)
        say("cannot build that lead: %s", exc)
        say("")
        say("RESULT: the lead could not be used")
        write_report(args.report, out)
        return 2

    say("=" * 72)
    say("%s", describe_set(leads))
    if args.show:
        for lead in leads:
            say("   -- %s", lead["label"])
            for key in sorted(PARAM_TYPES):
                if key in lead:
                    say("      %-16s %s", key, lead[key])
            if "x" in lead:
                say("      %-16s %.2f", "x (absolute)", lead["x"])
        say("")
        say("RESULT: OK")
        write_report(args.report, out)
        return 0

    meshes = load_meshes(args.stl_dir)
    corridor = json.load(open(args.corridor))
    foramen = (
        load_foramen_corridors(args.foramen) if os.path.exists(args.foramen) else None
    )

    rc = 0
    try:
        outcome = validate_set(
            leads,
            meshes=meshes,
            corridor=corridor,
            foramen=foramen,
            max_samples=args.max_samples,
            stl_dir=args.stl_dir,
        )
    except ValueError as exc:
        say("   could not validate: %s", exc)
        say("")
        say("RESULT: the lead could not be measured")
        write_report(args.report, out)
        return 2

    say("")
    for line in report_set_lines(args.tag, leads, outcome)[2:]:
        say("%s", line)
    if not outcome["ok"]:
        rc = 1

    if args.compare_rado:
        for lead, result in zip(leads, outcome["leads"]):
            if lead["type"] != "drg" or not is_built(result):
                continue
            say("")
            say("-" * 72)
            for line in compare_rado(lead, result, args.stl_dir):
                say("%s", line)

    if args.export:
        if not outcome["ok"] and not args.force:
            say("")
            say(
                "NOT EXPORTED: this lead does not fit. Fix it, or pass --force if "
                "you know what you are doing."
            )
        else:
            outdir, files = export_set(args.tag, leads, outcome["leads"], args.out)
            say("")
            say("wrote %d files to %s", len(files), os.path.relpath(outdir, REPO))
    say("")

    say(
        "RESULT: %s",
        {0: "OK", 1: "the lead does not fit", 2: "the lead could not be used"}[rc],
    )
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
