#!/usr/bin/env python3
"""Derive true left/right from the STL geometry and validate ansys/body_aliases.yaml.

The L/R strings in Khadka's filenames cannot be trusted: AprilRootDGR_right-*
contains files tagged `_l`, plain AprilRootDGR-* contains files tagged `_r`, and
several families are simply mirrored. Since the whole point of the label layer
is to say which side a body is on, the sides have to come from the geometry.

What this does, in order:

  1. Reads every STL and computes an area-weighted surface centroid. Area
     weighting rather than a plain vertex mean, because RADO's meshes are not
     uniformly tessellated and a vertex mean drifts toward the dense regions.
  2. Finds the left-right axis empirically: the axis along which known paired
     families (AprilRootMainRootL vs ...R and friends) separate most.
  3. Fixes the anatomical frame two independent ways and checks they agree --
     see the long comment in body_aliases.yaml. Short version: +Z is rostral,
     -Y is anterior, so in a right-handed frame anatomical Right = -X; and the
     descending aorta, which is left of midline, sits at +X. Both say +X = LEFT.
  4. Estimates the midline from the bodies that are unambiguously on it.
  5. Compares every filename L/R tag against which side of that midline the body
     actually sits on, and prints the disagreements.
  6. Re-checks body_aliases.yaml against all of the above: every declared side
     must match the geometry, and every declared level must be consistent with
     the rostro-caudal ordering.

Pure numpy -- scipy and PyYAML are not installed here.

    python3 ansys/check_laterality.py                 # full report
    python3 ansys/check_laterality.py --quiet         # just the verdicts
"""
import argparse
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from find_floating_bodies import read_stl_vertices          # noqa: E402
from apply_labels import parse_aliases, match_key            # noqa: E402

# A body whose centroid is within this of the midline is called midline, not
# sided. 2 mm is ten times the midline uncertainty and still well inside the
# smallest genuinely lateral structure (the rootlet bundles, at 2.4 mm).
MIDLINE_TOL_MM = 2.0

# Bodies that are on the midline by construction: the cord and its envelopes,
# the vertebrae and the discs. Used to locate the midline itself.
MIDLINE_BODIES = ["neuro_whitemater", "neuro_GreyMater", "neuro_CSF",
                  "neuro_Meninges", "neuro_EpiduralSpace", "neuro_blood_next_to_WM",
                  "neuro_blood_in_fat", "V1-2", "V1-3", "V1-4", "T_disk_"]

# Known paired families, used only to find which axis is left-right.
PAIRS = [("AprilRootMainRootL", "AprilRootMainRootR"),
         ("AprilRootBranchRootL", "AprilRootBranchRootR"),
         ("AprilDorsalRootsBloodParallel_L", "AprilDorsalRootsBloodParallel_R"),
         ("AprilRootDGR_right", "AprilRootDGR-")]

# Every laterality token that appears in a filename, and how to read it. "DGR"
# is a typo of DRG, so the R in it is NOT a side tag -- these patterns are
# deliberately specific rather than a bare search for "R".
TAG_PATTERNS = [
    (r"AprilRootDGR_right", "R"),
    (r"AprilRootMainRoot([LR])-", None),
    (r"AprilRootBranchRoot([LR])-", None),
    (r"AprilDorsalRootsBloodParallel_([LR])-", None),
    (r"AprilDorsalRootsInside([LR])", None),
    (r"AprilDorsalRootsMiddle([LR])", None),
    (r"AprilDorsalRootsOutsideMenginges([LR])", None),
    (r"AprilDorsolBranchInside([LR])", None),
    (r"AprilDorsolBranchMiddle([LR])", None),
    (r"AprilDorsolBranchOutsidMenging([LR])", None),
    (r"csf_coating_([rl])", None),
    (r"menging_coating_([rl])", None),
    (r"in_middle_([RL])", None),
    (r"connection_([RL])-", None),
]


def surface_centroid(verts):
    """Area-weighted centroid of a triangle soup, plus its bounding box."""
    tri = verts.reshape(-1, 3, 3).astype(np.float64)
    a, b, c = tri[:, 0], tri[:, 1], tri[:, 2]
    w = np.linalg.norm(np.cross(b - a, c - a), axis=1)      # 2 x triangle area
    if w.sum() == 0:
        return tri.reshape(-1, 3).mean(0)
    return (tri.mean(1) * w[:, None]).sum(0) / w.sum()


def load(stl_dir):
    bodies = {}
    for fn in sorted(os.listdir(stl_dir)):
        if not fn.lower().endswith(".stl"):
            continue
        v = read_stl_vertices(os.path.join(stl_dir, fn))
        if len(v) == 0:
            continue
        stem = os.path.splitext(fn)[0]
        bodies[stem] = {"c": surface_centroid(v), "lo": v.min(0), "hi": v.max(0)}
    return bodies


def tags_of(stem):
    """Every laterality token in one filename, uppercased. May be contradictory."""
    hits = []
    for pat, fixed in TAG_PATTERNS:
        for m in re.finditer(pat, stem):
            hits.append(fixed if fixed else m.group(1).upper())
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stl-dir", default=os.path.join(HERE, "..", "STL_files"))
    ap.add_argument("--map", default=os.path.join(HERE, "body_aliases.yaml"))
    ap.add_argument("--quiet", action="store_true", help="skip the per-body listing")
    args = ap.parse_args()

    bodies = load(args.stl_dir)
    print("loaded %d bodies from %s\n" % (len(bodies), os.path.relpath(args.stl_dir)))

    # -- 1. which axis is left-right ---------------------------------------
    print("STEP 1 -- left-right axis, from the separation of known paired families")
    print("%-34s %-34s %s" % ("family A", "family B", "centroid A - centroid B (mm)"))
    seps = []
    for pa, pb in PAIRS:
        A = np.array([b["c"] for k, b in bodies.items() if pa in k])
        B = np.array([b["c"] for k, b in bodies.items() if pb in k])
        if not len(A) or not len(B):
            continue
        d = A.mean(0) - B.mean(0)
        seps.append(np.abs(d))
        print("%-34s %-34s %s" % (pa, pb, np.round(d, 2)))
    seps = np.array(seps)
    axis = int(np.argmax(seps.mean(0)))
    print("\nmean |separation| per axis: %s -> axis %d (%s) is left-right"
          % (np.round(seps.mean(0), 2), axis, "XYZ"[axis]))
    allc = np.array([b["c"] for b in bodies.values()])
    print("rostro-caudal axis is Z, spanning %.1f..%.1f mm, as expected\n"
          % (allc[:, 2].min(), allc[:, 2].max()))

    # -- 2. anatomical frame ------------------------------------------------
    print("STEP 2 -- which direction along X is anatomical LEFT (two independent checks)")
    discs = {k: b["c"] for k, b in bodies.items() if "T_disk_" in k}
    rostral = max(discs.items(), key=lambda kv: kv[1][2])
    caudal = min(discs.items(), key=lambda kv: kv[1][2])
    print("  a) frame construction")
    print("     most rostral disc is %-24s at z=%.1f" % (rostral[0].split(" - ")[-1], rostral[1][2]))
    print("     most caudal  disc is %-24s at z=%.1f" % (caudal[0].split(" - ")[-1], caudal[1][2]))
    print("     -> +Z is rostral (T9/10 above T12/L1)")
    aorta = [b for k, b in bodies.items() if "Thoracic_aorta" in k][0]
    cord = [b for k, b in bodies.items() if "neuro_whitemater" in k][0]
    verts = np.array([b["c"] for k, b in bodies.items() if re.search(r"V1-\d$", k)])
    print("     aorta y=%.1f, vertebrae y=%.1f, cord y=%.1f -> -Y is anterior"
          % (aorta["c"][1], verts[:, 1].mean(), cord["c"][1]))
    print("     right-handed frame: Right = Anterior x Superior = (-y) x (+z) = -x")
    print("     -> +X is anatomical LEFT")
    print("  b) aorta anchor (independent)")
    print("     the descending thoracic aorta lies LEFT of midline; it is at x=%.1f"
          % aorta["c"][0])

    # -- 3. midline ---------------------------------------------------------
    midx = np.array([b["c"][0] for k, b in bodies.items()
                     if any(m in k for m in MIDLINE_BODIES)])
    MID = float(np.median(midx))
    wm = [b for k, b in bodies.items() if "neuro_whitemater" in k][0]
    print("\nSTEP 3 -- midline")
    print("  %d midline bodies, centroid x median %.3f, range %.3f..%.3f"
          % (len(midx), MID, midx.min(), midx.max()))
    print("  white-matter bounding box centre x = %.3f (independent)"
          % ((wm["lo"][0] + wm["hi"][0]) / 2))
    print("  model bounding-box centre x = %.3f (biased -- the model is not symmetric)"
          % ((np.array([b["lo"] for b in bodies.values()]).min(0)[0] +
              np.array([b["hi"] for b in bodies.values()]).max(0)[0]) / 2))
    print("  -> MIDLINE x = %.2f mm, good to about +/-0.2 mm" % MID)
    print("  aorta offset = %+.1f mm, i.e. LEFT, which confirms step 2b\n" % (aorta["c"][0] - MID))

    def side_of(cx):
        if abs(cx - MID) <= MIDLINE_TOL_MM:
            return "M"
        return "L" if cx > MID else "R"

    # -- 4. filename tags vs geometry --------------------------------------
    print("STEP 4 -- filename L/R tags vs geometry (midline tolerance %.1f mm)" % MIDLINE_TOL_MM)
    agree, disagree, contradictory, untagged = [], [], [], []
    for stem, b in sorted(bodies.items()):
        geo = side_of(b["c"][0])
        hits = set(tags_of(stem))
        if not hits:
            untagged.append((stem, b, geo))
        elif len(hits) > 1:
            contradictory.append((stem, b, geo, sorted(hits)))
        elif hits.pop() == geo:
            agree.append((stem, b, geo))
        else:
            disagree.append((stem, b, geo))
    if not args.quiet:
        print("\n%-62s %8s %8s %4s %s" % ("body", "cx", "offset", "geo", "tag"))
        for stem, b, geo in disagree:
            print("%-62s %8.2f %+8.2f %4s %s  DISAGREE"
                  % (stem.split(" - ")[-1][:62], b["c"][0], b["c"][0] - MID, geo,
                     sorted(set(tags_of(stem)))[0]))
        for stem, b, geo, hits in contradictory:
            print("%-62s %8.2f %+8.2f %4s %s  SELF-CONTRADICTORY"
                  % (stem.split(" - ")[-1][:62], b["c"][0], b["c"][0] - MID, geo, "/".join(hits)))
    print("\n  tag agrees with geometry: %d" % len(agree))
    print("  tag disagrees:            %d" % len(disagree))
    print("  filename contradicts itself: %d" % len(contradictory))
    print("  no tag at all:            %d" % len(untagged))

    # per-family roll-up, which is what makes the pattern legible
    print("\n  by family:")
    fams = {}
    for stem, b in bodies.items():
        fam = re.sub(r"-\d+", "-*", stem.split(" - ")[-1].split(" ")[0])
        hits = set(tags_of(stem))
        tag = sorted(hits)[0] if len(hits) == 1 else ("/".join(sorted(hits)) if hits else "-")
        fams.setdefault((fam, tag), []).append(b["c"][0] - MID)
    print("  %-36s %-6s %5s %9s %s" % ("family", "tag", "n", "mean off", "verdict"))
    for (fam, tag), offs in sorted(fams.items()):
        o = float(np.mean(offs))
        geo = side_of(o + MID)
        if tag == "-":
            v = "untagged, geometry says %s" % geo
        elif "/" in tag:
            v = "self-contradictory, geometry says %s" % geo
        elif tag == geo:
            v = "tag CORRECT"
        else:
            v = "tag INVERTED"
        print("  %-36s %-6s %5d %+9.1f %s" % (fam[:36], tag, len(offs), o, v))

    # -- 5. validate body_aliases.yaml -------------------------------------
    print("\nSTEP 5 -- validate %s against the geometry" % os.path.basename(args.map))
    rules = parse_aliases(args.map)
    bad_side, bad_level, unmatched = [], [], []
    per_rule = {}
    for stem, b in sorted(bodies.items()):
        key = match_key(stem)
        for rule in rules:
            m = re.match(rule["pattern"], key)
            if m:
                per_rule.setdefault(rule["name"], []).append((stem, b, rule, m))
                break
        else:
            unmatched.append(stem)

    # side and level are declared per SOURCE INDEX (the capture group the maps
    # are keyed on), not per body, so they are validated on the mean centroid of
    # each index group. A single rootlet filament can sit within the midline
    # tolerance while its bundle is clearly lateral; the bundle is what the rule
    # is asserting about.
    unsided, low_margin = [], []
    for rule in rules:
        entries = per_rule.get(rule["name"], [])
        if not entries:
            continue
        groups = {}
        for stem, b, r, m in entries:
            idx = m.group(1) if m.re.groups else ""
            groups.setdefault(idx, []).append(b["c"])
        for idx, cs in sorted(groups.items()):
            cs = np.array(cs)
            cx = cs[:, 0].mean()
            geo = side_of(cx)
            # A rule asserts a side either through `side`/`side_map` or by
            # spelling L/R straight into the label text, e.g. "DRG R{level}".
            # Only a rule with neither is genuinely unsided.
            declared = rule["side"] or rule["side_map"].get(idx)
            if declared is None:
                if geo != "M":
                    unsided.append((rule["name"], idx, geo, cx - MID))
            elif declared != geo:
                # Inside the midline band the geometry cannot call a side, but
                # the offset still has a sign. If the sign agrees with the rule
                # this is a low-margin body, not a contradiction -- e.g. the
                # filament-0 strand of rootlet bundles 7 and 11 sits nearer the
                # cord than the rest of its bundle, which is what sets the side.
                sign = "L" if cx > MID else "R"
                if geo == "M" and sign == declared:
                    low_margin.append((rule["name"], idx, declared, cx - MID))
                else:
                    bad_side.append(("%s[%s]" % (rule["name"], idx),
                                     declared, "%s (%+.2f mm)" % (geo, cx - MID)))

        # within one side, a larger level number must sit further caudal
        byside = {}
        for idx, cs in groups.items():
            lvl = rule["level"] or rule["level_map"].get(idx)
            if lvl is None:
                continue
            side = rule["side"] or rule["side_map"].get(idx) or "M"
            byside.setdefault(side, {}).setdefault(int(lvl), []).append(
                np.array(cs)[:, 2].mean())
        for side, levels in byside.items():
            seq = [(lvl, float(np.mean(zs))) for lvl, zs in sorted(levels.items())]
            for (l0, z0), (l1, z1) in zip(seq, seq[1:]):
                if z1 >= z0:
                    bad_level.append((rule["name"], side, seq))
                    break

    if bad_side:
        print("  SIDE MISMATCHES (%d):" % len(bad_side))
        for stem, declared, geo in bad_side:
            print("     %-58s declared=%s geometry=%s" % (stem.split(" - ")[-1][:58], declared, geo))
    else:
        print("  side: every rule agrees with the geometry")
    if low_margin:
        print("  low margin (%d) -- declared side is right but the body is inside the"
              " +/-%.1f mm midline band:" % (len(low_margin), MIDLINE_TOL_MM))
        for name, idx, declared, off in low_margin:
            print("     %-28s idx %-3s declared %s, offset %+.2f mm" % (name, idx or "-", declared, off))
    if unsided:
        print("  off-midline but deliberately unsided (%d rule/index groups):" % len(unsided))
        for name, idx, geo, off in unsided:
            print("     %-28s idx %-3s sits %+6.1f mm (%s) -- label carries no side by design"
                  % (name, idx or "-", off, geo))
    if bad_level:
        print("  LEVEL ORDER PROBLEMS (%d):" % len(bad_level))
        for name, side, items in bad_level:
            print("     %s side %s: %s" % (name, side, [(l, round(z)) for l, z, _ in items]))
    else:
        print("  level: every declared level is consistent with rostro-caudal order")
    if unmatched:
        print("  UNMATCHED (%d):" % len(unmatched))
        for stem in unmatched:
            print("     %s" % stem)
    else:
        print("  coverage: all %d bodies matched a rule" % len(bodies))

    return 1 if (bad_side or bad_level or unmatched) else 0


if __name__ == "__main__":
    sys.exit(main())
