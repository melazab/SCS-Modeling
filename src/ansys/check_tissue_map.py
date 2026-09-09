#!/usr/bin/env python3
"""Validate src/ansys/tissue_map.yaml against the actual STL filenames.

Every body must map to exactly one tissue. An unmatched body would silently get
a default material in Ansys, and a body matching two patterns means the map is
ambiguous. Both are worth catching before a multi-hour solve, not after.

Parsed without PyYAML so the script runs on a bare interpreter -- the file is a
small, regular subset of YAML, so it is parsed directly.

    python3 check_tissue_map.py [--stl-dir STL_files] [--map src/ansys/tissue_map.yaml]
"""
import argparse
import fnmatch
import os
import re
import sys


def parse_tissue_map(path):
    """Return [{name, ed_material, color_rgb, patterns}, ...] in file order.

    color_rgb is read for src/freecad/apply_colors.py, which imports this
    parser so the two scripts can never disagree about which body is which
    tissue. It is [r, g, b] 0-255, or None when the entry has no colour.
    """
    tissues, cur = [], None
    in_patterns = False
    for raw in open(path):
        line = raw.rstrip("\n")
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        m = re.match(r"^\s*-\s+name:\s*(\S+)", line)
        if m:
            cur = {"name": m.group(1), "ed_material": None,
                   "color_rgb": None, "patterns": []}
            tissues.append(cur)
            in_patterns = False
            continue
        if cur is None:
            continue
        m = re.match(r'^\s*ed_material:\s*"?([^"]*)"?\s*$', line)
        if m:
            cur["ed_material"] = m.group(1).strip()
            in_patterns = False
            continue
        m = re.match(r"^\s*color_rgb:\s*\[(.*)\]\s*$", line)
        if m:
            cur["color_rgb"] = [int(v) for v in m.group(1).split(",")]
            in_patterns = False
            continue
        m = re.match(r"^\s*patterns:\s*\[(.*)\]\s*$", line)
        if m:                                   # inline list
            body = m.group(1).strip()
            if body:
                cur["patterns"] = [p.strip().strip('"').strip("'")
                                   for p in body.split(",") if p.strip()]
            in_patterns = False
            continue
        if re.match(r"^\s*patterns:\s*$", line):
            in_patterns = True
            continue
        if in_patterns:
            # NB: list entries may carry a trailing "# comment" -- strip it, or
            # the entry is silently dropped (this bit me: "connection*" has one).
            m = re.match(r'^\s*-\s*(.+?)\s*$', line)
            if m:
                val = m.group(1)
                val = re.sub(r'\s+#.*$', '', val).strip().strip('"').strip("'")
                if val:
                    cur["patterns"].append(val)
                continue
            in_patterns = False
    return tissues


def matches(stem, pattern):
    """Case-insensitive match of one map pattern against a filename stem.

    Every RADO filename is prefixed with "T8-10 - ", so a glob has to be matched
    as a SUBSTRING, not against the whole string -- otherwise "T_disk_*" never
    matches "T8-10 - T_disk_10_11x-1".
    """
    p, s = pattern.lower(), stem.lower()
    if "*" in p or "?" in p:
        if not p.startswith("*"):
            p = "*" + p
        if not p.endswith("*"):
            p = p + "*"
        return fnmatch.fnmatch(s, p)
    return p in s


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument("--stl-dir", default=os.path.join(here, "..", "..", "STL_files"))
    ap.add_argument("--map", default=os.path.join(here, "tissue_map.yaml"))
    args = ap.parse_args()

    tissues = parse_tissue_map(args.map)
    print("parsed %d tissue entries from %s" % (len(tissues), os.path.relpath(args.map)))

    names = sorted(f for f in os.listdir(args.stl_dir)
                   if f.lower().endswith(".stl"))
    print("found %d STL files in %s\n" % (len(names), os.path.relpath(args.stl_dir)))

    counts = {t["name"]: 0 for t in tissues}
    unmatched, ambiguous = [], []

    for fn in names:
        stem = os.path.splitext(fn)[0]
        # patterns are matched case-insensitively against the filename stem,
        # as substrings unless they contain an explicit glob character
        hits = []
        for t in tissues:
            for pat in t["patterns"]:
                if matches(stem, pat):
                    hits.append(t["name"])
                    break
        if not hits:
            unmatched.append(stem)
        else:
            counts[hits[0]] += 1          # first match wins, per the map's rule
            if len(hits) > 1:
                ambiguous.append((stem, hits))

    print("%-24s %-24s %s" % ("TISSUE", "ENG-DATA MATERIAL", "BODIES"))
    total = 0
    for t in tissues:
        n = counts[t["name"]]
        total += n
        flag = "" if n else "   <-- no bodies matched"
        print("%-24s %-24s %5d%s" % (t["name"], t["ed_material"] or "-", n, flag))
    print("%-24s %-24s %5d" % ("TOTAL", "", total))

    if ambiguous:
        print("\nAMBIGUOUS (%d) -- matched more than one tissue; first wins:" % len(ambiguous))
        for stem, hits in ambiguous[:40]:
            print("   %s -> %s" % (stem, hits))

    if unmatched:
        print("\nUNMATCHED (%d) -- these would get a default material:" % len(unmatched))
        for stem in unmatched:
            print("   %s" % stem)
        return 1
    print("\nAll %d bodies mapped." % len(names))
    return 0


if __name__ == "__main__":
    sys.exit(main())
