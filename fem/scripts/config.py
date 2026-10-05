"""Shared configuration for the FreeCAD-side (gmsh + P1 FEM) SCS pipeline."""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
STL = os.path.join(ROOT, "STL_files")
LEAD = os.path.join(ROOT, "fem", "leads", "fem_dorsal_T10")
OUT = os.path.join(ROOT, "fem", "out")
LEAD_RUNS = os.path.join(OUT, "lead_runs")

TISSUE_MAP_YAML = os.path.join(ROOT, "src", "ansys", "tissue_map.yaml")

# The FIVE PRIMARY COMPARTMENTS, one STL each. This is NOT the classification
# set any more (see TISSUE_BODIES below, which covers all ~240 anatomical
# STLs) -- it is the short list of single, well-defined shells that two other
# things still need a ONE path for, and only for those two things:
#
#   * build_mesh.build()'s size-field clouds ("dura"/"white"/"canal" are
#     graded against exactly these surfaces) and mesh_cost.CostModel, which
#     mirrors them. Grading against 240 bodies would be a different -- and
#     far more expensive -- size field, and mesh_cost.py's RAW_TO_ACTUAL
#     calibration was measured against THIS one.
#   * scs_montage.py's tissue_surface(), which samples the solved field on
#     the dura/white/grey surfaces for the Potential Visualizer's 3D view.
#
# Every entry here also appears in TISSUE_BODIES under the same key (as the
# first, largest body of that class), so nothing is classified twice or
# defined twice.
TISSUE_STL = {
    "epidural": os.path.join(STL, "T8-10 - neuro_EpiduralSpace-1.STL"),
    "dura":     os.path.join(STL, "T8-10 - neuro_Meninges-1.STL"),
    "csf":      os.path.join(STL, "T8-10 - neuro_CSF-1.STL"),
    "white":    os.path.join(STL, "T8-10 - neuro_whitemater-1.STL"),
    "grey":     os.path.join(STL, "T8-10 - neuro_GreyMater-1.STL"),
}
CONTACT_STL = {i: os.path.join(LEAD, f"SCS Lead Electrode {i}.stl") for i in range(1, 9)}
INSULATOR_STL = os.path.join(LEAD, "SCS Lead Insulator.stl")


# ---------------------------------------------------------------------------
# THE FULL TISSUE SET (2026-09-15) -- every anatomical STL, not five
# compartments
# ---------------------------------------------------------------------------
# Until today this pipeline classified tets against exactly five shells
# (epidural, dura, CSF, white, grey) plus the lead, and DELETED every tet that
# fell outside all of them. fem/README.md's "What is approximate" item 3
# documented that as a deliberate simplification. Mohamed's instruction:
# "For now, just mesh everything. Put it as a TODO for later."
#
# NOTHING HERE NAMES A TISSUE, A CONDUCTIVITY OR A COLOUR OF ITS OWN. The
# filename->tissue-class mapping, the conductivities and the RGB colours all
# come from src/ansys/tissue_map.yaml at import time, through that directory's
# own parser (src/ansys/check_tissue_map.parse_tissue_map, the same one
# src/freecad/apply_colors.py and SCS_Mesh_Generator.FCMacro already use). The
# ONLY thing this file adds is the PRIORITY ORDER and the short pipeline names.
#
# WHY SHORT NAMES AT ALL: "white"/"grey"/"csf"/"dura"/"epidural" are baked into
# solution.npz's `order` array, analyze.py, plot_slices.py, scs_montage.py,
# crosscheck_elmer.py, dura_leak_report() and both FreeCAD macros. Renaming
# them to tissue_map's own spellings would churn all of that for nothing, so
# the five keep their names and the new classes get short names in the same
# style. ORDER_TO_TISSUE_CLASS is the one place the two vocabularies meet.
#
# ---------------------------------------------------------------------------
# PRIORITY, AND WHY IT IS THIS ORDER (first match wins in classify_tets)
# ---------------------------------------------------------------------------
# With five well-nested compartments ordering barely mattered. It matters now:
# RADO's bodies genuinely overlap, because they are a CAD assembly in which a
# containing body is often NOT hollowed out around what sits inside it. Counts
# below are TETS OF THE 66 939-TET COARSE MESH claimed by both classes,
# measured 2026-09-15 (a full all-pairs membership pass, not an assumption):
#
#     vertebra  ^ epidural  4856     vertebra ^ csf      2436
#     dura      ^ epidural  3339     dura     ^ csf      3080
#     dura      ^ vertebra  1655     root     ^ dura      566
#     root      ^ csf        484     blood    ^ epidural  365
#     disc      ^ epidural   289     blood    ^ csf       250
#     root      ^ white      250     disc     ^ vertebra  243
#     root      ^ epidural   174     root     ^ grey      157
#     blood     ^ vertebra   128     blood    ^ dura      118
#     ... plus 10 more pairs under 30 tets each.
#
# 15 555 tets (23 %) are inside exactly two classes and 1 078 inside three or
# more, so the order below decides the tissue of roughly a quarter of the mesh.
# The rule applied, in one sentence: THE MOST SPECIFIC BODY WINS WHERE IT IS A
# GENUINE INCLUSION, AND THE ENCLOSING BULK WINS WHERE THE SMALLER BODY MERELY
# BLEEDS INTO IT. Concretely, in priority order:
#
#  1. grey, white -- the cord. Nothing carves the cord: it is the best-defined
#     geometry in the set and the one the published numbers are measured on.
#     grey before white because the horns sit inside the columns (2 tets
#     overlap, unchanged from the original five-compartment order).
#  2. drg, root, blood, sympathetic -- discrete small structures threaded
#     THROUGH the fluid/fat spaces. A rootlet inside the dural sac is nerve,
#     not CSF (root ^ csf 484); a radicular artery in epidural fat is vessel,
#     not fat (blood ^ epidural 365); a root piercing the dura wall is root,
#     not dura (root ^ dura 566). drg before root because the ganglion's
#     "football_in_middle" core sits inside the root bundle (6 tets).
#     RADO's root bundles come in Inside / Middle / OutsidMenging layers:
#     tissue_map puts Inside and Middle in nerve_root and OutsidMenging in
#     meninges_dura (its *Menging* pattern wins), so "root beats dura" is
#     exactly what makes the layering come out right -- the core is nerve and
#     the sheath around it is dura. See the TODO below about Middle.
#  3. csf, dura, epidural -- the bulk fluid/membrane/fat layers of the canal,
#     innermost first, exactly as the original five-compartment order had
#     them. csf before dura matters more than it used to: the root sheaths
#     (class dura) run through the subarachnoid space, and without this a
#     sheath would carve dura out of the CSF (dura ^ csf 3080).
#  4. disc, vertebra -- the bony envelope, LAST. This is not cosmetic:
#     RADO's three vertebra STLs are SOLID bodies whose spinal canal is not
#     carved out, so they contain 17 926 of the coarse mesh's 66 939 tets,
#     including 4 856 of the epidural space and 2 436 of the CSF. Anywhere but
#     last, "vertebra" would swallow the cord. disc before vertebra for the
#     243 tets they share.
#
# ---------------------------------------------------------------------------
# TODO (UNRESOLVED, AFFECTS THE FIELD RIGHT NEXT TO THE ELECTRODES)
# ---------------------------------------------------------------------------
# tissue_map.yaml itself flags an open question it has not settled: RADO models
# each root as Inside / Middle / OutsidMenging, mirroring the DRG's
# in_middle / csf_coating / menging_coating triple. If that analogy holds, the
# "...Middle..." bodies are CSF (1.7 S/m), not nerve (0.1432 S/m) -- a 12x
# conductivity difference in bodies that sit right beside the contacts.
# tissue_map.yaml currently assigns them nerve_root, and THIS FILE FOLLOWS
# tissue_map.yaml, deliberately: nothing here hardcodes an answer to a question
# the source of truth has not answered. Settle it in tissue_map.yaml (by
# measurement or by asking whoever built the CAD), not here. Also recorded in
# fem/README.md "What is approximate".
# ---------------------------------------------------------------------------

# pipeline name -> tissue_map.yaml `name`. KEY ORDER IS THE PRIORITY ORDER.
ANATOMY_CLASS = [
    ("grey",        "grey_matter"),
    ("white",       "white_matter"),
    ("drg",         "dorsal_root_ganglion"),
    ("root",        "nerve_root"),
    ("blood",       "blood_vessel"),
    ("sympathetic", "sympathetic_chain"),
    ("csf",         "csf"),
    ("dura",        "meninges_dura"),
    ("epidural",    "epidural_space"),
    ("disc",        "intervertebral_disc"),
    ("vertebra",    "vertebra"),
]

# tissue_map classes deliberately NOT in ANATOMY_CLASS:
#   electrode_contact / lead_insulation -- LEAD HARDWARE, not anatomy. The
#     "SCS Lead Electrode 1-4.stl"/"SCS Lead Insulator.stl" files sitting in
#     STL_files/ are RADO's own 4-contact DRG lead. This pipeline's lead is
#     whichever one the caller passes in (the frozen fem/leads/fem_dorsal_T10,
#     or the live tree lead live_lead.py tessellates), and it already gets
#     "contactN"/"insulator" entries of its own at the head of order_for().
#     Classifying RADO's spare lead as tissue too would silently drop a SECOND
#     electrode array -- four floating 1e4 S/m metal bodies -- into every
#     solve. TODO: model RADO's own DRG lead when a DRG montage is actually
#     wanted; it should come in as a lead, through contact_stl=, not as tissue.
#   soft_tissue -- has no `patterns` and matches no STL (the thorax envelope
#     is not in RADO's STL export). Nothing to classify.


def _tissue_sigma(path):
    """{tissue_map name: sigma_S_per_m} straight out of the yaml.

    check_tissue_map.parse_tissue_map() reads name/ed_material/color_rgb/
    patterns but not the conductivity, and that file belongs to the Ansys
    side (read-only from here), so this adds the one missing field by the
    same line-oriented approach on the same file. It does NOT restate any
    value: the numbers live in tissue_map.yaml and nowhere else.
    """
    out, cur = {}, None
    with open(path) as fh:
        for line in fh:
            m = re.match(r"^\s*-\s+name:\s*(\S+)", line)
            if m:
                cur = m.group(1)
                continue
            m = re.match(r"^\s*sigma_S_per_m:\s*([-+0-9.eE]+)\s*$", line)
            if m and cur is not None:
                out[cur] = float(m.group(1))
    return out


def _load_tissue_map():
    """(entries, sigma_by_class) from src/ansys/tissue_map.yaml, using that
    directory's own parser so there is exactly one filename->tissue matcher
    in the repo."""
    ansys_dir = os.path.join(ROOT, "src", "ansys")
    if ansys_dir not in sys.path:
        sys.path.insert(0, ansys_dir)
    from check_tissue_map import parse_tissue_map
    return parse_tissue_map(TISSUE_MAP_YAML), _tissue_sigma(TISSUE_MAP_YAML)


def _classify_stl_files():
    """{pipeline name: [STL path, ...]} for every file in STL_files/ that
    tissue_map.yaml matches to one of ANATOMY_CLASS's classes, first pattern
    match wins exactly as check_tissue_map.py's own main() does. Files
    matching a class NOT in ANATOMY_CLASS (the lead hardware) are skipped;
    a file matching nothing at all is returned separately so a caller can
    say so instead of silently losing it."""
    entries, _ = _load_tissue_map()
    from check_tissue_map import matches
    want = dict(ANATOMY_CLASS)
    by_class = {k: [] for k, _ in ANATOMY_CLASS}
    unmatched = []
    for fn in sorted(os.listdir(STL)):
        if not fn.lower().endswith(".stl"):
            continue
        stem = os.path.splitext(fn)[0]
        hit = None
        for t in entries:
            if any(matches(stem, p) for p in t["patterns"]):
                hit = t["name"]
                break
        if hit is None:
            unmatched.append(fn)
            continue
        for name, cls in ANATOMY_CLASS:
            if cls == hit:
                by_class[name].append(os.path.join(STL, fn))
                break
    return by_class, unmatched


_BY_CLASS, UNMATCHED_STL = _classify_stl_files()

# {pipeline name: [STL path, ...]} -- the CLASSIFICATION set. A class with no
# body in STL_files/ is dropped rather than left as an entry that can never be
# hit (it would still cost a name slot in every `order` and every report).
TISSUE_BODIES = {k: v for k, v in _BY_CLASS.items() if v}

# Priority-ordered anatomy names actually present, i.e. ANATOMY_CLASS filtered
# to what STL_files/ really contains.
ANATOMY_ORDER = tuple(k for k, _ in ANATOMY_CLASS if _BY_CLASS[k])

# pipeline name -> tissue_map.yaml class name, for anything that needs the
# colour or the Engineering Data material (SCS_Mesh_Generator.FCMacro's
# per-facet tissue colours read this, then look the RGB up in the yaml).
ORDER_TO_TISSUE_CLASS = {k: c for k, c in ANATOMY_CLASS if _BY_CLASS[k]}


def order_for(contact_stl):
    """Tissue/contact classification order for an ARBITRARY lead: one
    "contactN" entry per key actually present in `contact_stl` (ascending),
    then insulator, then the anatomy in ANATOMY_ORDER's priority order --
    generalizes ORDER below (which is just order_for(CONTACT_STL), the
    frozen 8-contact lead) to any contact count N. Every function that
    classifies tets or interprets a `lab`/`sigma` array by tissue NAME
    (assign_and_solve.classify_tets/sigma_of/dura_leak_report/
    solve_basis_fields) takes an explicit `order` parameter that defaults to
    this, computed from whichever contact_stl dict was actually passed in --
    so a 4-contact live lead gets a 4-long contact block here, not a fixed
    8-long one with 4 entries pointing at STLs that don't exist. First hit
    wins; metal before insulator (contacts are recessed in the insulator
    body), then the anatomy.

    The anatomy block is ANATOMY_ORDER: eleven classes covering all ~240
    anatomical STLs. `lab` values are POSITIONS INTO THIS LIST, so a stored
    solution.npz or mesh_report.json must be read against the `order` array
    saved beside it, never against today's config.ORDER."""
    return (["contact%d" % i for i in sorted(contact_stl)]
            + ["insulator"] + list(ANATOMY_ORDER))


# The frozen field-of-record lead's own order (8 contacts) -- kept as a
# module-level constant for every EXISTING caller that reads config.ORDER
# directly assuming the frozen lead (main()'s own default args, this
# module's own CONTACT_STL/INSULATOR_STL). Equal to order_for(CONTACT_STL);
# NOT the order used for an arbitrary live lead -- see order_for() above.
ORDER = order_for(CONTACT_STL)


def _build_sigma():
    """Electrical conductivity, S/m, per pipeline tissue name.

    Values are read from src/ansys/tissue_map.yaml's sigma_S_per_m and never
    restated here, so the two cannot drift.

    Two entries are not anatomy classes:
      insulator  -- the lead's own insulation (tissue_map lead_insulation).
      background -- the fill used only by the SCS_BACKGROUND=1 variant run
                    for tets outside every body; tissue_map's vertebra value,
                    as it always was.
    """
    _, sig = _load_tissue_map()
    out = {name: sig[cls] for name, cls in ORDER_TO_TISSUE_CLASS.items()}
    out["insulator"] = sig["lead_insulation"]
    out["background"] = sig["vertebra"]
    return out


SIGMA = _build_sigma()
# tissue_map gives metal electrode 4.0e6 S/m. A 1e8 conductivity ratio against
# epidural fat wrecks the conditioning of the linear system for no physical
# gain: at 1e4 S/m a contact is already equipotential to ~1e-5 of the driving
# voltage. SIGMA_METAL_TRUE is kept so the clamp can be audited.
SIGMA_METAL_TRUE = 4.0e6
SIGMA_METAL = 1.0e4

# Bipolar drive, as specified: source contact, sink contact, current in amperes.
SOURCE_CONTACT = 3
SINK_CONTACT = 5
CURRENT_A = 1.0


# ---------------------------------------------------------------------------
# Lead-directory resolution, shared by SCS_Mesh_Generator.FCMacro and
# SCS_Potential_Visualizer.FCMacro. Pure path/filesystem logic, no FreeCAD
# dependency, so it's usable from fem/.venv subprocesses too.
#
# There is no "pick which lead" concept: the lead is whatever sits in the live
# document tree. live_lead.py finds it by object name and never reads an export
# directory. Do not reintroduce a lead-selection dropdown here.
#
# lead_stls() and is_single_lead_dir() are generic over "a directory with N
# contact STLs plus one insulator STL" and do not care where it came from;
# live_lead.py hands them fem/out/lead_runs/_live/, the scratch directory it
# tessellates the live lead into.
# ---------------------------------------------------------------------------
def detect_contact_count(lead_dir):
    """How many contacts `lead_dir`'s flat layout actually has, by counting
    "SCS Lead Electrode <i>.stl" files present for i = 1, 2, 3, ... stopping
    at the first gap. This is what lets a reader that doesn't already know
    N (mesh_preview.py's LEAD_DIR, solve_lead.py's --lead-dir -- both come
    from a CLI/JSON argument, not from live_lead.py's own in-memory count)
    discover it directly from what's on disk, instead of assuming 8. A
    live_lead.py caller that DOES already know N (it just tessellated and
    is about to write exactly N contact STLs) can skip this probe by
    passing n_contacts= to lead_stls() directly."""
    n = 0
    while os.path.isfile(os.path.join(lead_dir, "SCS Lead Electrode %d.stl" % (n + 1))):
        n += 1
    return n


def lead_stls(lead_dir=None, n_contacts=None):
    """(contact_stl dict, insulator_stl path) for `lead_dir`'s single-lead
    flat layout, or for the frozen field-of-record LEAD (always 8 contacts)
    if lead_dir is None. For an explicit lead_dir, `n_contacts` contacts are
    assumed present (1..n_contacts); if omitted it is discovered from the
    directory itself via detect_contact_count(). Never mutates
    CONTACT_STL/INSULATOR_STL -- callers thread the returned dict/path
    through explicitly (build_mesh.build(), assign_and_solve.
    classify_tets()/solve_basis_fields()) rather than pointing this module's
    own globals somewhere else out from under any other code that imported
    them at module-load time."""
    if lead_dir is None:
        d, n = LEAD, 8
    else:
        d = lead_dir
        n = detect_contact_count(d) if n_contacts is None else n_contacts
    contact = {i: os.path.join(d, "SCS Lead Electrode %d.stl" % i) for i in range(1, n + 1)}
    insulator = os.path.join(d, "SCS Lead Insulator.stl")
    return contact, insulator


def is_single_lead_dir(path):
    """True if `path` directly holds one lead's N contact STLs (N >= 1) +
    insulator (the flat layout build_lead_config.export_set() writes for a
    SINGLE lead, and the layout fem/leads/fem_dorsal_T10 has always had --
    that one happens to have N=8, but nothing here requires it). A
    MULTI-lead export instead holds per-lead subdirectories (lead1_dorsal/,
    ...) plus a leads.txt manifest and no flat STLs at top level -- that
    layout is OUT OF SCOPE for the interactive mesh/solve pipeline (see
    SCS_Mesh_Generator.FCMacro's module docstring); this returns False for
    it so callers can detect and reject it with a clear message instead of
    silently meshing zero contacts. Also False for a directory with an
    insulator but NO contact STLs at all -- an empty `contact` dict would
    otherwise vacuously satisfy all(... for p in contact.values())."""
    contact, insulator = lead_stls(path)
    return bool(contact) and os.path.isfile(insulator) and all(os.path.isfile(p) for p in contact.values())
