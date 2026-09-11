#!/bin/bash
# ONE VARIABLE: replace the 878 electrode constraint equations with coupled
# DOF sets. Everything else -- materials, pinball, the strip, the +/-1 A
# bipolar injection -- is identical to job 3802339, so any change in the
# pivots is attributable to CE vs CP and nothing else.
#
# WHY. Job 3802339 put the maximum pivot, 7.733301104E+16, at node 3149 --
# which is an electrode CE master, carrying 505 of the deck's 878 `volt`
# constraint equations. Stripping a third of the model moved that number not
# at all, and the ECC sweep moved it not at all across four orders of
# magnitude. The conditioning lives in the electrode representation.
#
# Every one of those 878 lines is the same trivial statement:
#
#     ce,next,0.,3149,volt,-1.,3150,volt,1.     "node 3150 has 3149's voltage"
#
# That is what CP,,VOLT means, and the deck contains 878 `ce,` lines and ZERO
# `cp,` lines. The difference is not cosmetic:
#
#   CE  a general linear constraint. MAPDL may carry it with a Lagrange
#       multiplier, which adds a row and makes the assembled system a
#       saddle-point problem -- indefinite, and negative pivots are a normal
#       feature of those. Ours fails on exactly that: a negative pivot on a
#       VOLT degree of freedom.
#   CP  a coupled set. Implemented by ELIMINATION: the coupled DOFs collapse
#       into a single unknown. No multiplier, no extra row, and the matrix
#       gets SMALLER rather than harder.
#
# The deck has exactly two electrodes, nodes 786 (373 slaves) and 3149 (505
# slaves), and already injects the bipolar current as a boundary condition:
#
#     f,3149,amps,-1000000000.      -1 A   (uMKS: picoamps)
#     f,786,amps,1000000000.        +1 A
#
# So current injection needs no change at all -- it is already a surface
# condition, not a property of the meshed metal. Only the equipotential
# representation changes here. Removing the meshed platinum is the NEXT step
# and is deliberately not combined with this one; two changes at once would
# make an ambiguous result.
#
# CPNGEN rather than CP: CP takes at most 17 nodes per command and redefines
# the set when the same NSET is repeated, whereas CPNGEN ADDS a node to a
# set, which is what a 506-node electrode needs.
#
#SBATCH --job-name=scs_cp
#SBATCH --account=tlv
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=400G
#SBATCH --time=04:00:00
#SBATCH --output=slurm_%j.out
set -uo pipefail
module load ansys/25.1
export OMP_NUM_THREADS=1
export KMP_AFFINITY=disabled

NRANK=12
PINB=${PINB:-500}
SRC=$HOME/Documents/gatech/scs_model_bipolar_current.dat
WORK=${PFSDIR:-${TMPDIR:-/tmp}}/scs_cp_$SLURM_JOB_ID
DEST=$HOME/Documents/SCS-Modeling/ansys_training/bigdeck_cp/out_$SLURM_JOB_ID
mkdir -p "$WORK" "$DEST"
echo "=== $(date) on $(hostname); CE->CP; stripped; PINB=-${PINB} ==="
free -g | head -2

echo "--- patching deck ---"
time sed -e 's/\r$//' \
    -e "s/^rmod,cid,6,0\./rmod,cid,6,-${PINB}./" \
    -e "s/^rmod,tid,6,0\./rmod,tid,6,-${PINB}./" \
    -e 's|^/fclean|!/fclean  ! disabled: keep the .rst|' \
    "$SRC" > "$WORK/scs_model.dat"

ok=1
grep -q "^rmod,cid,6,-${PINB}\."  "$WORK/scs_model.dat" || { echo "FAIL cid pinball"; ok=0; }
grep -q "^MP,RSVX,1,2\.5e-13,"    "$WORK/scs_model.dat" || { echo "FAIL metal changed!"; ok=0; }
grep -q "^MP,RSVX,5,0\.05,"       "$WORK/scs_model.dat" || { echo "FAIL insulation changed!"; ok=0; }
[ "$ok" -eq 1 ] || { echo "=== aborting: patch incomplete ==="; exit 1; }
echo "materials confirmed at original Khadka values"

python3 - "$WORK/scs_model.dat" <<'PY'
import io, re, sys
p = sys.argv[1]
src = io.open(p, encoding="latin-1").read().split("\n")

# Gather the equipotential sets: master -> [slaves]. Every ce line in this deck
# is the two-term form "0 = -1*master + 1*slave"; anything else is left alone
# and counted, so a deck that is not this shape cannot be silently mangled.
pat = re.compile(r"^ce,next,0\.,(\d+),volt,-1\.,(\d+),volt,1\.\s*$", re.I)
sets, other, kept = {}, 0, []
for ln in src:
    m = pat.match(ln)
    if m:
        sets.setdefault(m.group(1), []).append(m.group(2))
    else:
        if ln.lower().startswith("ce,"):
            other += 1
        kept.append(ln)

if other:
    sys.exit("ABORT: %d ce lines are not the simple equipotential form" % other)
if not sets:
    sys.exit("ABORT: no equipotential ce lines found")

lines = ["", "!  ---- CE -> CP: equipotential electrodes as coupled DOF sets ----"]
for i, (master, slaves) in enumerate(sorted(sets.items(), key=lambda kv: int(kv[0])), 1):
    lines.append("/com, ==== CP set %d: master %s, %d slaves ====" % (i, master, len(slaves)))
    lines.append("cpngen,%d,volt,%s" % (i, master))
    for n in slaves:
        lines.append("cpngen,%d,volt,%s" % (i, n))
lines.append("cpintf,volt,1e-6   ! nothing to couple by proximity; harmless, documents intent")
lines.append("!  ---- end CE -> CP ----")

# Insert where the ce block was: immediately before the first /solu that follows
# the deck's own boundary conditions. Anchor on the current injection, which
# must stay exactly as it is.
anchor = None
for i, ln in enumerate(kept):
    if re.match(r"^f,\d+,amps,", ln, re.I):
        anchor = i
        break
if anchor is None:
    sys.exit("ABORT: no f,,amps line found to anchor on")
out = kept[:anchor] + lines + kept[anchor:]
io.open(p, "w", encoding="latin-1").write("\n".join(out))
print("converted %d electrodes, %d nodes, to coupled sets"
      % (len(sets), sum(len(v) + 1 for v in sets.values())))
PY

echo "--- verify conversion ---"
ok=1
grep -q "^cpngen,1,volt,"            "$WORK/scs_model.dat" || { echo "FAIL no cpngen"; ok=0; }
grep -qi "^ce,next,0\.,[0-9]*,volt"  "$WORK/scs_model.dat" && { echo "FAIL ce lines remain"; ok=0; }
grep -qi "^f,3149,amps,-1000000000\."  "$WORK/scs_model.dat" || { echo "FAIL current injection changed"; ok=0; }
grep -qi "^f,786,amps,1000000000\."    "$WORK/scs_model.dat" || { echo "FAIL current injection changed"; ok=0; }
[ "$ok" -eq 1 ] || { echo "=== aborting: conversion incomplete ==="; exit 1; }
echo -n "cpngen lines: "; grep -c "^cpngen," "$WORK/scs_model.dat"
echo "bipolar +/-1 A injection unchanged"

# same strip as job 3802339, so this differs from it in exactly one respect
python3 - "$WORK/scs_model.dat" <<'PY'
import io, sys
p = sys.argv[1]
strip = """
/prep7
allsel,all
*get,_e_before,elem,0,count
/com, ==== STRIP before: elements %_e_before% ====
cmsel,s,VERTEBRA
cmsel,a,VASCULATURE
cmsel,a,INTERVERTEBRAL_DISK
nsle,s
esln,s,0
edele,all
allsel,all
esel,all
nsle,s
nsel,inve
ndele,all
allsel,all
*get,_e_after,elem,0,count
/com, ==== STRIP after: elements %_e_after% ====
finish
/solu
"""
src = io.open(p, encoding="latin-1").read()
i = src.find("\nsolve")
if i < 0:
    sys.exit("ABORT: no solve found")
io.open(p, "w", encoding="latin-1").write(src[:i] + "\n" + strip + src[i:])
print("strip injected")
PY
grep -q "STRIP after" "$WORK/scs_model.dat" || { echo "=== aborting: strip missing ==="; exit 1; }

cd "$WORK"
echo "=== solve start: $(date) ==="
time ansys251 -b -dis -mpi intelmpi -np "$NRANK" -i scs_model.dat -o scs_model.out
echo "=== exit $? at $(date) ==="
if grep -qaE "SOLVE is not a recognized|No results file available" scs_model.out; then
    echo "FAIL: the solve was skipped."
fi
echo "--- markers ---"; grep -aE "^ ==== " scs_model.out | head
echo "--- was CP eliminated or multiplied? ---"
grep -aiE "coupled|constraint equation|lagrange|eliminat" scs_model.out | head -10
echo "--- pivots ---"
grep -aiE "maximum pivot|minimum pivot|\*\*\* error|SOLUTION IS DONE" scs_model.out | tail -10
cp -v scs_model.out file0.err "$DEST"/ 2>/dev/null
