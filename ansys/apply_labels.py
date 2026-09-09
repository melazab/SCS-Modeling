#!/usr/bin/env python3
"""Set a readable Label on every object in NBF_RADO-SCS.FCStd from body_aliases.yaml.

Khadka's STL filenames are the provenance link to the published RADO-SCS model
and are referenced by the Ansys .wbpj / .mechdb / .pmdb files and the 4.5 GB
MAPDL deck, so they are never renamed. A FreeCAD object has an immutable Name
(auto-derived from the filename, what everything keys on) and a free-text Label
(what the model tree shows, nothing depends on it). This script only ever writes
Label, so it cannot break the Ansys chain: reverting it changes nothing but the
strings in the tree.

    python3 ansys/apply_labels.py --dry-run                     # no FreeCAD needed
    freecadcmd ansys/apply_labels.py                            # writes the document
    APPLY_LABELS_ARGS="--clean-core" freecadcmd ansys/apply_labels.py

freecadcmd owns the command line: it consumes every flag itself, and passing
flags after its --pass makes it skip the script entirely. So under freecadcmd
the script's own options come from the APPLY_LABELS_ARGS environment variable.
Bare `freecadcmd ansys/apply_labels.py` applies every default, which is the
normal way to run it.

Three things about freecadcmd that shape this script:

  - freecadcmd IMPORTS the script rather than exec'ing it, so __name__ is not
    "__main__" and a plain main-guard would leave the script doing nothing at
    all, silently. See run_as_freecadcmd_script() at the bottom.

  - print() output is swallowed, so everything is written to a report file
    (--report, default ansys/apply_labels_report.txt) which the caller cats.
  - freecadcmd segfaults on exit AFTER the script has finished and the document
    is safely saved. That is harmless but it means the exit code is meaningless:
    verify from the report file, which ends with a "RESULT:" line. The segfault
    also drops a `core` file in the cwd. It lands after this process is gone, so
    the script can only clear the PREVIOUS run's core (--clean-core, done first
    thing); the caller still has to rm the one this run leaves behind.

The script is idempotent. It matches on Name (via the STL filename it was
derived from), never on the current Label, so re-running after a successful run
is a no-op and re-running after a partial run finishes the job.
"""
import argparse
import os
import re
import shlex
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
DEFAULT_DOC = os.path.join(REPO, "NBF_RADO-SCS.FCStd")
DEFAULT_STL = os.path.join(REPO, "STL_files")
DEFAULT_MAP = os.path.join(HERE, "body_aliases.yaml")
DEFAULT_REPORT = os.path.join(HERE, "apply_labels_report.txt")

PREFIX = "T8-10 - "          # stripped before matching; every RADO body has it


# --------------------------------------------------------------------------
# body_aliases.yaml reader
#
# PyYAML is not installed here, and the file is a small regular subset of YAML,
# so it is parsed directly -- same approach as check_tissue_map.py. Only the
# keys this script uses are recognised; anything else in a rule is ignored.
# --------------------------------------------------------------------------
def parse_aliases(path):
    """Return [{name, pattern, label, side, side_map, level, level_map}, ...] in file order."""
    rules, cur = [], None
    for raw in open(path, encoding="utf-8"):
        line = raw.rstrip("\n")
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        m = re.match(r"^\s*-\s+name:\s*(\S+)", line)
        if m:
            cur = {"name": m.group(1), "pattern": None, "label": None,
                   "side": None, "side_map": {}, "level": None, "level_map": {}}
            rules.append(cur)
            continue
        if cur is None:
            continue
        # scalar values are quoted with ' or " in this file; strip whichever
        m = re.match(r"^\s*(pattern|label|side|level):\s*(.+?)\s*$", line)
        if m:
            key, val = m.group(1), m.group(2)
            if len(val) >= 2 and val[0] == val[-1] and val[0] in "'\"":
                val = val[1:-1]
            cur[key] = val
            continue
        m = re.match(r"^\s*(side_map|level_map):\s*(.+?)\s*$", line)
        if m:
            key, val = m.group(1), m.group(2).strip("'\"")
            pairs = {}
            for item in val.split(","):
                if ":" in item:
                    k, v = item.split(":", 1)
                    pairs[k.strip()] = v.strip()
            cur[key] = pairs
            continue
    bad = [r["name"] for r in rules if not r["pattern"] or not r["label"]]
    if bad:
        raise ValueError("rules missing pattern or label: %s" % ", ".join(bad))
    return rules


def script_args():
    """Arguments meant for this script, under either interpreter.

    Run as `python3 apply_labels.py --dry-run`, sys.argv is the usual thing. Run
    as `freecadcmd apply_labels.py`, sys.argv[0] is "freecadcmd" and FreeCAD has
    swallowed the whole command line -- and `--pass`, which looks like the way
    to forward flags, actually makes freecadcmd skip running the script at all.
    So options come from APPLY_LABELS_ARGS instead.
    """
    if os.path.basename(sys.argv[0]).startswith("freecad"):
        return shlex.split(os.environ.get("APPLY_LABELS_ARGS", ""))
    return sys.argv[1:]


def match_key(stem):
    """Filename stem as the rules see it: the shared 'T8-10 - ' prefix removed."""
    return stem[len(PREFIX):] if stem.startswith(PREFIX) else stem


def render(rule, m):
    """Build a label from a rule and its regex match, or raise if a lookup is missing."""
    out = rule["label"]
    idx = m.group(1) if m.re.groups >= 1 else None

    if "{side}" in out:
        side = rule["side"] or rule["side_map"].get(idx)
        if side is None:
            raise KeyError("rule %s has no side for index %r" % (rule["name"], idx))
        out = out.replace("{side}", side)

    if "{level}" in out:
        level = rule["level"] or rule["level_map"].get(idx)
        if level is None:
            raise KeyError("rule %s has no level for index %r" % (rule["name"], idx))
        out = out.replace("{level}", str(level))

    for g in range(1, m.re.groups + 1):
        out = out.replace("{%d}" % g, m.group(g) or "")
    return out


def resolve(stem, rules):
    """Return (label, rule_name) for one filename stem, or (None, None). First match wins."""
    key = match_key(stem)
    for rule in rules:
        m = re.match(rule["pattern"], key)
        if m:
            return render(rule, m), rule["name"]
    return None, None


def sanitize(stem):
    """Reproduce the object Name FreeCAD derives from an STL filename.

    Mesh import replaces every character that is not [A-Za-z0-9_] with '_'. The
    mapping is 1:1 over this STL set (checked below), which is what lets the
    script key on Name and stay idempotent once the Labels have changed.
    """
    return re.sub(r"[^A-Za-z0-9_]", "_", stem)


def build_index(stl_dir, rules):
    """Return (name -> (stem, label, rule)), plus the list of stems that matched nothing."""
    stems = sorted(os.path.splitext(f)[0] for f in os.listdir(stl_dir)
                   if f.lower().endswith(".stl"))
    index, unmatched, collisions = {}, [], []
    for stem in stems:
        label, rule = resolve(stem, rules)
        if label is None:
            unmatched.append(stem)
            continue
        name = sanitize(stem)
        if name in index:
            collisions.append((name, index[name][0], stem))
        index[name] = (stem, label, rule)
    return index, unmatched, collisions, stems


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--doc", default=DEFAULT_DOC)
    ap.add_argument("--stl-dir", default=DEFAULT_STL)
    ap.add_argument("--map", default=DEFAULT_MAP)
    ap.add_argument("--report", default=DEFAULT_REPORT)
    ap.add_argument("--dry-run", action="store_true",
                    help="report what would change; touches neither the document nor FreeCAD")
    ap.add_argument("--no-backup", action="store_true",
                    help="skip the .FCStd copy (only for a re-run over an already-backed-up doc)")
    ap.add_argument("--clean-core", action="store_true",
                    help="delete the `core` file the previous freecadcmd run's exit segfault left")
    args = ap.parse_args(script_args())

    out = []
    def say(fmt, *a):
        out.append(fmt % a if a else fmt)

    if args.clean_core and os.path.isfile("core"):
        os.remove("core")
        say("removed the core file left by the previous freecadcmd exit segfault")

    rules = parse_aliases(args.map)
    index, unmatched, collisions, stems = build_index(args.stl_dir, rules)
    say("mode:      %s", "DRY RUN" if args.dry_run else "WRITE")
    say("rules:     %d from %s", len(rules), os.path.relpath(args.map, REPO))
    say("STL files: %d in %s", len(stems), os.path.relpath(args.stl_dir, REPO))
    say("document:  %s", args.doc)
    say("")

    used = {}
    for _, (stem, label, rule) in index.items():
        used.setdefault(rule, []).append(label)
    say("%-26s %6s  %s", "RULE", "BODIES", "EXAMPLE LABEL")
    for rule in rules:
        hits = used.get(rule["name"], [])
        flag = "" if hits else "   <-- matched nothing"
        say("%-26s %6d  %s%s", rule["name"], len(hits),
            sorted(hits)[0] if hits else "-", flag)
    say("%-26s %6d", "TOTAL", sum(len(v) for v in used.values()))
    say("")

    dup = {}
    for _, (stem, label, _r) in index.items():
        dup.setdefault(label, []).append(stem)
    clashes = {k: v for k, v in dup.items() if len(v) > 1}
    if clashes:
        say("DUPLICATE LABELS (%d) -- readable but ambiguous in the tree:", len(clashes))
        for label, group in sorted(clashes.items()):
            say("   %-46s %s", label, group)
        say("")
    if collisions:
        say("NAME COLLISIONS (%d) -- two STLs sanitize to one object Name:", len(collisions))
        for name, a, b in collisions:
            say("   %s  <-  %s | %s", name, a, b)
        say("")
    if unmatched:
        say("UNMATCHED (%d) -- no rule produced a label:", len(unmatched))
        for stem in unmatched:
            say("   %s", stem)
        say("")

    if args.dry_run:
        # Without FreeCAD there is no "before" to diff against, so show the
        # mapping itself. The filename stem IS the label FreeCAD starts with.
        say("WOULD SET (stem -> label), first 25 of %d:", len(index))
        for name in sorted(index)[:25]:
            stem, label, _r = index[name]
            say("   %-62s -> %s", stem[:62], label)
        say("")
        say("RESULT: dry run OK, %d bodies would be labelled, %d unmatched",
            len(index), len(unmatched))
        write_report(args.report, out)
        return 0 if not unmatched else 1

    # -- write path: FreeCAD only imported here so --dry-run works anywhere ---
    import shutil
    import FreeCAD

    if not args.no_backup:
        bak = args.doc + ".prelabel.bak"
        n = 1
        while os.path.exists(bak):
            n += 1
            bak = "%s.prelabel%d.bak" % (args.doc, n)
        shutil.copy2(args.doc, bak)
        say("backup:    %s", os.path.basename(bak))

    doc = FreeCAD.openDocument(args.doc)
    say("opened:    %s, %d objects", doc.Name, len(doc.Objects))
    say("")

    changed, already, missing_rule, samples = 0, 0, [], []
    for obj in doc.Objects:
        entry = index.get(obj.Name)
        if entry is None:
            missing_rule.append((obj.Name, obj.Label))
            continue
        stem, label, _r = entry
        before = obj.Label
        if before == label:
            already += 1
            continue
        obj.Label = label
        changed += 1
        samples.append((obj.Name, before, label))

    say("changed:   %d", changed)
    say("unchanged: %d (already correct)", already)
    if missing_rule:
        say("")
        say("OBJECTS WITH NO RULE (%d) -- Label left alone:", len(missing_rule))
        for name, label in missing_rule:
            say("   %-58s Label=%r", name[:58], label)
    say("")
    say("BEFORE -> AFTER (first 25 changes):")
    for name, before, after in samples[:25]:
        say("   %-62s", before[:62])
        say("       -> %s", after)

    if changed:
        doc.save()
        say("")
        say("saved:     %s", args.doc)
    else:
        say("")
        say("no changes, document not rewritten")
    FreeCAD.closeDocument(doc.Name)

    say("")
    say("RESULT: %d changed, %d already correct, %d objects without a rule, "
        "%d STL stems unmatched", changed, already, len(missing_rule), len(unmatched))
    write_report(args.report, out)
    return 0


def write_report(path, lines):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def run_as_freecadcmd_script():
    """True when freecadcmd was handed THIS file to run.

    freecadcmd does not exec a script, it imports it as a module, so __name__ is
    "apply_labels" and the usual __main__ guard never fires -- the script loads,
    defines everything, and silently does nothing. Comparing the path freecadcmd
    was given against __file__ is what makes `freecadcmd ansys/apply_labels.py`
    actually run, without also firing when check_laterality.py imports this
    module for its parser.
    """
    if len(sys.argv) < 2 or not os.path.basename(sys.argv[0]).startswith("freecad"):
        return False
    return os.path.abspath(sys.argv[1]) == os.path.abspath(__file__)


if __name__ == "__main__":
    sys.exit(main())
elif run_as_freecadcmd_script():
    # No sys.exit here: raising SystemExit out of an imported module makes
    # freecadcmd print a traceback over a run that actually succeeded.
    main()
