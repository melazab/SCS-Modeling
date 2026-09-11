#!/bin/bash
# DEBUG RUN: full model minus vertebrae, vasculature and intervertebral discs.
#
# Purpose is diagnosis, not physics. Job 3800027 showed the pivot RATIO
# (5.55e15 / 4.6e-06) is what now blocks the solve, and that ECC is not the
# lever -- three values four orders apart gave bit-identical pivots. Removing
# a third of the elements will not by itself fix conditioning, but it makes
# every subsequent experiment cheap: the full model costs ~40 min and 420 GB
# per attempt, which is far too slow a loop for the number of things left to
# try.
#
# NOT A PHYSICS MODEL. Voiding conductive bodies changes the current
# distribution:
#   vertebra          0.04 S/m   nearly insulating -- mild approximation
#   intervertebral    0.6  S/m   conductive, in the path
#   vasculature       0.66 S/m   conductive, threads past the roots
# Discs and vasculature leave insulating cavities where conductive tissue was,
# so nothing from this run should be quoted as a result. It exists to make the
# matrix smaller and the turnaround faster.
#
# SELECT BY COMPONENT, NOT BY MATERIAL. The deck defines named components
# (cm,VERTEBRA / cm,VASCULATURE / cm,INTERVERTEBRAL_DISK / cm,EPIDURAL_SPACE
# and others) and that is what makes this safe. Selecting by material number
# would be a trap: vertebra and epidural space BOTH carry material 6, because
# Engineering Data gives them the same 25 ohm*m, so `esel,s,mat,,6` would
# delete the epidural space -- the one compartment the lead sits in.
#
# CONTACT ELEMENTS MUST GO WITH THEM. Deleting a body's solid elements while
# leaving the CONTA174/TARGE170 pairs on its surface would leave those pairs
# referencing nodes with no stiffness behind them -- which is precisely the
# electrically floating condition that caused the original
# "insufficiently constrained model" failures. So the deletion is done as:
# select the bodies, take their nodes, select EVERY element touching those
# nodes (esln with key 0 catches the contact elements), delete all of it, then
# delete the nodes nothing references any more.
#
#SBATCH --job-name=scs_stripped
#SBATCH --account=tlv
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=600G
#SBATCH --time=06:00:00
#SBATCH --output=slurm_%j.out
set -uo pipefail
module load ansys/25.1
export OMP_NUM_THREADS=1
export KMP_AFFINITY=disabled

NRANK=12
PINB=${PINB:-500}
SRC=$HOME/Documents/gatech/scs_model_bipolar_current.dat
WORK=${PFSDIR:-${TMPDIR:-/tmp}}/scs_stripped_$SLURM_JOB_ID
DEST=$HOME/Documents/SCS-Modeling/ansys_training/bigdeck_stripped/out_$SLURM_JOB_ID
mkdir -p "$WORK" "$DEST"
echo "=== $(date) on $(hostname); stripped; PINB=-${PINB} um; ranks=$NRANK ==="
free -g | head -2

echo "--- patching deck (materials untouched) ---"
time sed -e 's/\r$//' \
    -e "s/^rmod,cid,6,0\./rmod,cid,6,-${PINB}./" \
    -e "s/^rmod,tid,6,0\./rmod,tid,6,-${PINB}./" \
    -e 's|^/fclean|!/fclean  ! disabled: keep the .rst|' \
    "$SRC" > "$WORK/scs_model.dat"

echo "--- verify patch ---"
ok=1
grep -q "^rmod,cid,6,-${PINB}\."  "$WORK/scs_model.dat" || { echo "FAIL cid pinball"; ok=0; }
grep -q "^!/fclean"               "$WORK/scs_model.dat" || { echo "FAIL fclean"; ok=0; }
grep -q "^MP,RSVX,1,2\.5e-13,"    "$WORK/scs_model.dat" || { echo "FAIL metal changed!"; ok=0; }
grep -q "^MP,RSVX,5,0\.05,"       "$WORK/scs_model.dat" || { echo "FAIL insulation changed!"; ok=0; }
for c in VERTEBRA VASCULATURE INTERVERTEBRAL_DISK EPIDURAL_SPACE; do
  grep -q "^cm,$c" "$WORK/scs_model.dat" || { echo "FAIL component $c missing"; ok=0; }
done
[ "$ok" -eq 1 ] || { echo "=== aborting: patch incomplete ==="; exit 1; }
echo "materials confirmed at original Khadka values; all components present"

# Inject the strip immediately before the first solve, inside /PREP7.
python3 - "$WORK/scs_model.dat" <<'PY'
import sys, io
p = sys.argv[1]
strip = """
!  ---- DEBUG STRIP: vertebrae, vasculature, intervertebral discs ----
/prep7
allsel,all
*get,_e_before,elem,0,count
*get,_n_before,node,0,count
/com, ==== STRIP before: elements %_e_before% nodes %_n_before% ====
cmsel,s,VERTEBRA
cmsel,a,VASCULATURE
cmsel,a,INTERVERTEBRAL_DISK
*get,_e_target,elem,0,count
/com, ==== STRIP targeting %_e_target% body elements ====
nsle,s
esln,s,0                 ! every element touching them, contact pairs included
*get,_e_kill,elem,0,count
/com, ==== STRIP deleting %_e_kill% elements (bodies + their contacts) ====
edele,all
allsel,all
esel,all
nsle,s                   ! nodes still referenced by a surviving element
nsel,inve                ! the orphans
*get,_n_kill,node,0,count
/com, ==== STRIP deleting %_n_kill% orphaned nodes ====
ndele,all
allsel,all
*get,_e_after,elem,0,count
*get,_n_after,node,0,count
/com, ==== STRIP after: elements %_e_after% nodes %_n_after% ====
finish
/solu
!  ---- END DEBUG STRIP ----
!  The /solu above is load-bearing. This block is injected immediately before
!  the deck's own `solve`, which sits inside a /SOLU section. Entering /prep7
!  to delete elements leaves that section, and a bare `finish` drops MAPDL to
!  BEGIN level -- where `solve` is not a valid command and is SILENTLY IGNORED
!  with only a warning. Job 3802282 did exactly that: the strip deleted
!  everything correctly, MAPDL skipped the solve, and the run "completed" in
!  5:48 with no results file and exit code 0.
"""
src = io.open(p, encoding="latin-1").read()
i = src.find("\nsolve")
if i < 0:
    sys.exit("no solve found - aborting")
out = src[:i] + "\n" + strip + src[i:]
io.open(p, "w", encoding="latin-1").write(out)
print("strip injected before the first solve")
PY
grep -q "END DEBUG STRIP" "$WORK/scs_model.dat" || { echo "=== aborting: strip not injected ==="; exit 1; }

cat >> "$WORK/scs_model.dat" <<'POST'

/post1
set,last
nsort,volt
*get,vmx,sort,0,max
*get,vmn,sort,0,min
/com, ==== STRIPPED VOLT range: min %vmn% max %vmx% ====
allsel,all
fini
POST

cd "$WORK"
echo "=== solve start: $(date) ==="
time ansys251 -b -dis -mpi intelmpi -np "$NRANK" -i scs_model.dat -o scs_model.out
echo "=== exit $? at $(date) ==="
echo "--- did the solve actually run? ---"
if grep -qaE "SOLVE is not a recognized|No results file available" scs_model.out; then
    echo "FAIL: the solve was skipped -- MAPDL was not in /SOLU when `solve` arrived."
fi
echo "--- strip markers ---"; grep -aE "^ ==== STRIP" scs_model.out
echo "--- volt ---";          grep -aE "^ ==== STRIPPED" scs_model.out
echo "--- pivots ---"
grep -aiE "maximum pivot|minimum pivot|pivot in absolute|\*\*\* error|SOLUTION IS DONE" scs_model.out | tail -20
cp -v scs_model.out file0.err "$DEST"/ 2>/dev/null
{ echo "work: $WORK"; ls -la "$WORK" | head -30; } > "$DEST/WORK_LOCATION.txt"
