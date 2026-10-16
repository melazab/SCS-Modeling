"""Which STL set the lead tools measure and check against.

The Lead Designer never reads anatomy from the open document: build_lead_config,
measure_corridor, measure_foramen and rado_drg_lead all load STL files from
disk by name. This module is the one place that says which directory.

    rado_up4 (default)  fem/out/geometry/rado_up4, the extended T6-T12 model as
                        separate bodies -- the same bodies NBF_RADO-SCS.FCStd
                        shows, and the source of the fused set the FEM meshes
                        (fem/scripts/config.py). It keeps the per-foramen
                        blocks and ganglia a DRG lead is checked against, which
                        the fused set merges.
    rado                STL_files/, RADO's own T10-T12 bodies.

It follows fem/scripts/config.py's switch, so the two cannot disagree:
SCS_GEOMETRY=rado selects RADO's own set for both, anything else the extended
one. The extended set is generated (src/extend/extend_rado.py), not committed.
"""
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))

SETS = {
    "rado_up4": os.path.join(REPO, "fem", "out", "geometry", "rado_up4"),
    "rado": os.path.join(REPO, "STL_files"),
}
NAME = "rado" if os.environ.get("SCS_GEOMETRY") == "rado" else "rado_up4"
STL_DIR = SETS[NAME]
PREFIX = "T8-10 - "


def require(stl_dir=None):
    stl_dir = stl_dir or STL_DIR
    if not os.path.isdir(stl_dir):
        raise IOError("anatomy STL set missing: %s -- build it with "
                      "fem/.venv/bin/python src/extend/extend_rado.py, or set "
                      "SCS_GEOMETRY=rado for RADO's own STLs" % stl_dir)
    return stl_dir


def stl_name(stem, stl_dir=None):
    """Filename of RADO body `stem` in this set: "T8-10 - <stem>.STL", or the
    extended model's single continuous body "T8-10 - <stem>__extended_upN.STL"
    that replaces it. Exactly one must exist."""
    stl_dir = require(stl_dir)
    pat = re.compile(r"^%s%s(__extended_up\d+)?\.stl$" % (re.escape(PREFIX), re.escape(stem)), re.I)
    hits = sorted(f for f in os.listdir(stl_dir) if pat.match(f))
    if len(hits) != 1:
        raise IOError("%d files for %r in %s" % (len(hits), stem, stl_dir))
    return hits[0]


def names_matching(regex, stl_dir=None):
    """Sorted filenames in the set whose stem (prefix removed) matches."""
    stl_dir = require(stl_dir)
    pat = re.compile(regex)
    out = []
    for f in sorted(os.listdir(stl_dir)):
        if not f.lower().endswith(".stl") or not f.startswith(PREFIX):
            continue
        if pat.search(os.path.splitext(f)[0][len(PREFIX):]):
            out.append(f)
    return out


def level_copy(stem):
    """(source stem, k) for a level copy "<stem>__up<k>", else (stem, 0)."""
    m = re.match(r"^(.*)__up(\d+)$", stem)
    return (m.group(1), int(m.group(2))) if m else (stem, 0)


# RADO's top vertebra is T10 in disc naming (lead_defaults.yaml, "Vertebral
# levels"), so copy k sits at level T(10 - k); the template labels the copied
# ganglia "DRG L1 -- ganglion core [T9 level]" and so on.
TOP_LEVEL = 10


def level_of_copy(k):
    return "T%d" % (TOP_LEVEL - k)
