#!/usr/bin/env python3
"""Rebuild the two committed leads and assert their STLs come back byte for byte.

WHY THIS FILE EXISTS
--------------------
src/freecad/generated_leads/dorsal_z110_8c/ and ventral_z110_8c/ hold 18 STLs
that are committed to the repository -- eight contacts and an insulator each.
They are the only lead geometry in this repo that anything downstream has ever
been pointed at, and they are the output of build_lead_config.py. So they are
also a fixture: if this module's sweep, its contact spacing, its insulator fuse
or its export ever quietly changes, rebuilding those leads produces different
bytes, and comparing the bytes is how you find out.

It has caught real drift twice. It is not a formality.

This used to be an entry in a `configs:` catalogue -- "export dorsal_z110_8c and
diff it". The catalogue is gone (see build_lead_config.py's docstring), so the
two parameter sets live HERE, written out in full, which is strictly better for
a test: they no longer depend on lead_defaults.yaml at all. Change a default
there -- the house diameter, the tail length -- and this test still builds the
same two leads and still passes, because it names every number itself. That is
the point of a fixture; a test that follows the defaults it is testing cannot
detect a change in them.

WHERE THE ODD z_center COMES FROM
---------------------------------
z_center is 110.43203353881836 and NOT the round 110.0 the directory name
suggests. Those STLs were built before any of this existed, by make_scs_lead.py
with no --z-center at all, and its default is the midpoint of the measured
corridor: 0.5 * (61.43203353881836 + 159.43203353881836). The "z110" in the
directory name is that number through "%.0f". Every digit is spelled out
deliberately -- at 110.0 the lead rebuilds 0.43 mm caudal of the committed STLs,
which is exactly the kind of silent divergence this test exists to catch.

AND A THIRD FIXTURE, WHICH IS NOT A BYTE COMPARISON
--------------------------------------------------
The DRG lead has no committed STLs to diff, and deliberately so: it is not
generated any more. It is RADO's own four-contact lead -- the five files in
STL_files/, which are already committed and already the source of truth --
rigidly transformed onto whichever ganglion is named (rado_drg_lead.py). So the
fixture it needs is not "do the bytes match a copy we kept", it is the stronger
claim the mechanism actually makes:

    1. ON RADO'S OWN GANGLION, WITH NO OFFSETS, THE TRANSFORM IS THE IDENTITY.
       Every vertex of all five bodies must come back BIT-IDENTICAL to the
       source STL. Not "within a tolerance" -- identical, because copying a
       lead onto the ganglion it is already on must not move it.

    2. ON ALL EIGHT TARGETS THE MAP IS RIGID AND THE MESH SURVIVES IT.
       Signed volume and surface area preserved exactly, zero open edges, and
       signed volume still POSITIVE -- that last one is the whole reason this
       case exists. Reaching a right-side ganglion from a left-side one is a
       reflection, which inverts every facet normal unless the vertex winding
       is reversed with it, and an inside-out STL breaks meshing downstream
       while looking perfectly fine in a viewport.

    3. THE DETERMINANT MATCHES THE SIDE. +1 within a side, -1 across, and the
       winding reversed exactly when it is -1.

    4. IT IS DETERMINISTIC. Placing twice gives byte-identical files.

None of that needs FreeCAD -- rado_drg_lead.py is numpy and struct -- so this
part runs even in a plain Python interpreter, and it runs FIRST, before the two
swept leads give up for want of a FreeCAD.

RUN IT
------
    freecadcmd src/freecad/test_regression.py     # all three fixtures
    python3 src/freecad/test_regression.py        # the DRG fixture only

Nothing is written into the repository: the rebuild goes to a temporary
directory which is removed afterwards, and the committed files are only read.
There is no document involved at any point, so none of the save/colour hazards
the rest of this directory is careful about apply here. Exit status is 0 if
everything that could be checked reproduced, 1 if anything did not.
"""
import filecmp
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
REPORT = os.path.join(HERE, "test_regression_report.txt")

# The defaults this test resolves against: EMPTY. Every number a lead needs is
# given below, so nothing is inherited and lead_defaults.yaml cannot move the
# fixture. `levels` is empty for the same reason -- these leads are positioned by
# an absolute z, not by a level name whose value could be re-derived.
CFG = {"defaults": {}, "epidural_defaults": {}, "drg_defaults": {},
       "levels": {}, "max_leads_per_type": 2}

Z_CENTER = 110.43203353881836

# The two fixtures: the parameter set, and the directory its STLs are committed
# in. An 8-contact 1.30 mm lead, 3 mm contacts 1 mm apart, 6 mm tails, on the
# midline -- one in the dorsal epidural space and one in the ventral.
CASES = [
    ("dorsal_z110_8c",
     {"type": "dorsal", "contacts": 8, "contact_length": 3.0, "gap": 1.0,
      "diameter": 1.3, "tail": 6.0, "x_offset": 0.0, "z_center": Z_CENTER}),
    ("ventral_z110_8c",
     {"type": "ventral", "contacts": 8, "contact_length": 3.0, "gap": 1.0,
      "diameter": 1.3, "tail": 6.0, "x_offset": 0.0, "z_center": Z_CENTER}),
]

# The nine files each case must reproduce. Named explicitly rather than globbed,
# so a case that exported six files and stopped fails loudly instead of passing
# on the six it managed.
EXPECTED = ["SCS Lead Electrode %d.stl" % i for i in range(1, 9)] \
    + ["SCS Lead Insulator.stl"]


def compare(case_dir, built_dir, say):
    """Byte-compare the nine STLs of one case. Returns True if all nine match."""
    ok = True
    for fn in EXPECTED:
        want = os.path.join(case_dir, fn)
        got = os.path.join(built_dir, fn)
        if not os.path.exists(want):
            say("   %-28s MISSING from the repository", fn)
            ok = False
            continue
        if not os.path.exists(got):
            say("   %-28s NOT BUILT", fn)
            ok = False
            continue
        # shallow=False: compare contents, not (size, mtime). A rebuild always
        # has a new mtime, so a shallow compare would report every file changed.
        same = filecmp.cmp(want, got, shallow=False)
        say("   %-28s %s  (%d bytes)", fn,
            "identical" if same else "DIFFERS", os.path.getsize(want))
        if not same:
            ok = False
    extra = sorted(set(os.listdir(built_dir)) - set(EXPECTED))
    if extra:
        say("   built %d file(s) nobody asked for: %s", len(extra), ", ".join(extra))
        ok = False
    return ok


def check_drg(say):
    """The DRG fixture. Returns True if it holds. Needs no FreeCAD.

    See the module docstring for what is being asserted and why it is not a
    byte comparison against a committed copy.
    """
    import numpy as np
    import rado_drg_lead as rdl

    ok = True
    say("=" * 72)
    say("DRG -- RADO's own lead, copied onto all eight ganglia")

    # 0. the constant is not stale: re-derive which ganglion RADO's lead is on.
    foramen = rdl.load_foramen()
    native, ranking = rdl.verify_native_target(foramen)
    say("   native ganglion re-derived from the meshes: %s (nearest core centroid "
        "%.2f mm; runner-up %s at %.2f mm)",
        native, ranking[0]["centroid_gap_mm"], ranking[1]["tag"],
        ranking[1]["centroid_gap_mm"])
    if native != rdl.NATIVE_TARGET:
        say("   FAIL: rado_drg_lead.NATIVE_TARGET says %s", rdl.NATIVE_TARGET)
        ok = False

    # 1. the identity case, bit for bit.
    source = {fn: tris for _k, _i, fn, tris in rdl.load_rado_bodies()}
    placed = rdl.place(rdl.NATIVE_TARGET, foramen)
    exact = True
    for _kind, _idx, fn, tris in placed["bodies"]:
        if tris.shape != source[fn].shape or not np.array_equal(tris, source[fn]):
            say("   FAIL: %s moved when copied onto its own ganglion", fn)
            exact = ok = False
    say("   %s onto %s with no offsets: all 5 bodies bit-identical to STL_files/",
        "PASS" if exact else "FAIL", rdl.NATIVE_TARGET)

    # 2 and 3. every target: rigid, watertight, outward, and the right parity.
    say("   %-4s %-6s %-9s %-9s %10s %10s %7s %s",
        "tgt", "side", "det", "winding", "vol", "area", "edges", "verdict")
    for tag in sorted(foramen["ganglia"]):
        p = rdl.place(tag, foramen)
        side = p["target_frame"]["side"]
        want_mirror = side != p["source_frame"]["side"]
        good = (p["orientation_ok"]
                and p["mirrored"] == want_mirror
                and abs(abs(p["det"]) - 1.0) < 1e-12
                and all(r["winding_reversed"] == want_mirror for r in p["orientation"]))
        worst = max(abs(abs(r["signed_volume_after"]) - abs(r["signed_volume_before"]))
                    for r in p["orientation"])
        worst_a = max(abs(r["area_after"] - r["area_before"]) for r in p["orientation"])
        edges = max(r["open_edges"] for r in p["orientation"])
        say("   %-4s %-6s %+9.6f %-9s %10.1e %10.1e %7d %s",
            tag, side, p["det"], "reversed" if p["mirrored"] else "kept",
            worst, worst_a, edges, "OK" if good else "FAILED")
        ok = ok and good

    # 4. determinism.
    a = rdl.place("R2", foramen)
    b = rdl.place("R2", foramen)
    same = all(np.array_equal(x[3], y[3]) for x, y in zip(a["bodies"], b["bodies"]))
    say("   placing R2 twice gives identical vertices: %s", "PASS" if same else "FAIL")
    ok = ok and same
    say("")
    return ok


def main():
    out = []

    def say(fmt, *a):
        out.append(fmt % a if a else fmt)

    sys.path.insert(0, HERE)
    import build_lead_config as blc

    drg_ok = True
    try:
        drg_ok = check_drg(say)
    except Exception as exc:                    # noqa: BLE001 -- report, do not crash
        say("DRG fixture could not run: %s: %s", type(exc).__name__, exc)
        say("")
        drg_ok = False

    say("rebuilding %d committed leads and comparing every byte", len(CASES))
    say("committed under %s", os.path.relpath(blc.DEFAULT_OUT, REPO))
    say("")

    try:
        blc._geometry_modules()
    except ImportError as exc:
        say("Building a SWEPT lead needs FreeCAD, and this interpreter has none "
            "(%s). The DRG fixture above does not, and its verdict stands.", exc)
        say("")
        say("    freecadcmd src/freecad/test_regression.py")
        say("")
        say("RESULT: %s", "DRG PASS, swept leads not checked" if drg_ok
            else "DRG FAIL")
        write(out)
        return 0 if drg_ok else 1

    meshes = blc.load_meshes()
    import json
    corridor = json.load(open(blc.DEFAULT_CORRIDOR))

    failed = []
    tmp = tempfile.mkdtemp(prefix="scs_regression_")
    try:
        for tag, params in CASES:
            case_dir = os.path.join(blc.DEFAULT_OUT, tag)
            say("=" * 72)
            say("%s", tag)

            lead = blc.resolve_lead(CFG, params)
            # The resolve must not have invented anything: every parameter the
            # lead carries is either one this file gave it or one derived from
            # those. A default leaking in here would mean the fixture is no
            # longer self-contained even though the bytes happened to match.
            derived = {"label", "index", "of", "x"}
            leaked = sorted(k for k in lead
                            if k not in params and k not in derived)
            if leaked:
                say("   FAIL: %s came from somewhere other than this file",
                    ", ".join(leaked))
                failed.append(tag)
                continue
            wrong = sorted(k for k, v in params.items() if lead[k] != v)
            if wrong:
                say("   FAIL: resolve changed %s", ", ".join(wrong))
                failed.append(tag)
                continue
            say("   %s", blc.describe(lead))

            outcome = blc.validate_set([lead], meshes=meshes, corridor=corridor)
            if not outcome["ok"]:
                # Not fatal to the byte comparison, but worth shouting about:
                # these two placements passed when they were committed.
                say("   NOTE: this lead no longer passes its own fit check")
            if not outcome["leads"][0].get("shapes"):
                say("   FAIL: the lead could not be built at all")
                failed.append(tag)
                continue

            built, _files = blc.export_set(tag, [lead], outcome["leads"],
                                           out_root=tmp)
            if not compare(case_dir, built, say):
                failed.append(tag)
            say("")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    say("=" * 72)
    if failed or not drg_ok:
        if not drg_ok:
            say("the DRG fixture FAILED -- see the top of this report")
        if failed:
            say("%d of %d swept leads did NOT reproduce: %s",
                len(failed), len(CASES), ", ".join(failed))
        say("")
        say("Something in the lead geometry has changed. Either the change is")
        say("wrong, or it is right and the committed STLs need regenerating --")
        say("but that is a decision, not a fixup, because everything downstream")
        say("of those files was meshed and solved from the old geometry.")
        say("")
        say("RESULT: FAIL")
        write(out)
        return 1
    say("all %d swept leads reproduced their %d committed STLs byte for byte, "
        "and the DRG fixture holds on all eight ganglia",
        len(CASES), len(CASES) * len(EXPECTED))
    say("")
    say("RESULT: PASS")
    write(out)
    return 0


def write(lines):
    with open(REPORT, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    sys.stdout.write("\n".join(lines) + "\n")


def run_as_freecadcmd_script():
    """True when a FreeCAD interpreter was handed THIS file; see apply_labels.py."""
    if len(sys.argv) < 2 or not os.path.basename(sys.argv[0]).lower().startswith("freecad"):
        return False
    return os.path.abspath(sys.argv[1]) == os.path.abspath(__file__)


if __name__ == "__main__":
    sys.exit(main())
elif run_as_freecadcmd_script():
    # No sys.exit: SystemExit out of an imported module makes freecadcmd print a
    # traceback over a run that succeeded.
    main()
