#!/bin/bash
# Full 20.8M-node model, ORIGINAL Khadka material values, single change:
# widen the contact pinball radius so the bonded electric contact pairs can
# actually close the gaps that exist in the geometry.
#
# WHY. The Sep 9 smp run (job 3796783, 331 GB, so memory was NOT the limit)
# failed with:
#
#     Distributed sparse solver maximum pivot = 1.92307552E+16 at node 7395321
#     Distributed sparse solver minimum pivot = -77.33 at node 7192526
#     minimum pivot in absolute value = 1.942856684E-07 at node 8711510 VOLT
#     *** ERROR *** There is at least 1 small equation solver pivot term ...
#                   Please check for an insufficiently constrained model.
#
# Two things that settles:
#   1. The max pivot is 1.92307552E+16 against 1.92280267E+16 in the ORIGINAL
#      August run -- identical to 5 significant figures despite an 11-order
#      material change. The metal/insulation contrast never drove the pivot
#      ratio. It comes from the contact conductance the deck sets itself:
#          *set,_maxCond,4.e+016
#          rmod,cid,19,_maxCond/_ASMDIAG      ! ECC
#      So the earlier "reconditioning" was both unacceptable and pointless. It
#      has been removed from every script; materials are Khadka's throughout.
#   2. "Insufficiently constrained model" with a ~1e-7 pivot is the signature of
#      an electrically FLOATING region.
#
# WHICH bodies float -- answered from the geometry, not guessed
# (src/freecad/find_floating_bodies.py builds a proximity graph over the 245 STL
# surfaces and finds its connected components):
#
#     gap tolerance      components   isolated bodies
#     0.05 mm            38           79
#     0.10 mm             9            8
#     0.20 mm             3            2   (Thoracic_aorta-1, connection_R-2)
#     0.30 mm             1            0   <-- fully connected
#
# So the model IS geometrically connected, but only once contact can bridge
# ~0.3 mm. The stragglers are the vascular chain -- thoracic aorta, the
# "connection" radicular vessels, the DorsalRootsBloodParallel segments -- which
# are separate STL bodies that merely abut. Blood is electrically continuous
# through them, so bonding them is the physically correct reading, not a hack.
#
# THE FIX. Every one of the 1116 contact pairs currently carries
#     rmod,cid,6,0.      ! PINB
#     rmod,tid,6,0.      ! PINB
# i.e. pinball left at the program default, evidently under the ~0.3 mm gap.
# Set it to an absolute 500 um instead (negative PINB = absolute distance; the
# deck is /units,uMKS so length is micrometres).
#
# Note this canNOT short unrelated structures: pinball only lets an EXISTING
# contact pair close a larger gap. Workbench created pairs from proximity, so
# there is no pair across e.g. the dura wall for a wider pinball to bridge.
#
# ECC SWEEP. Job 3797962 (pinball 500 um, direct sparse, smp) COMPLETED in
# 37 min on 418 GB and settled the previous question:
#
#     original 3796783        pinball 3797962
#     max pivot  1.92e16      max pivot  7.13e16
#     min pivot -77.33        min pivot -1.563
#     min |pivot| 1.94e-07    (gone)
#     "insufficiently         "large negative pivot" +
#      constrained model"     "extremely large pivot ratio"
#
# The 1.9e-07 tiny-pivot signature -- the electrically floating region -- is
# GONE. Widening the pinball to an absolute 500 um across all 2232 rmod lines
# did what it was meant to do. What remains is a different failure: a pivot
# RATIO of 7.13e16 / 1.56 ~ 4.6e16.
#
# And note the max pivot went UP, 1.92e16 -> 7.13e16. That is the tell. More
# contact pairs are now actually closing, so more of the deck's own contact
# conductance is being assembled into the matrix. The deck sets, once per pair:
#
#     *set,_maxCond,4.e+016
#     rmod,cid,19,_maxCond/_ASMDIAG      ! ECC
#
# 4e16 is essentially the max pivot, so the pivot IS the contact conductance.
#
# WHAT THIS IS NOT. ECC is not a material property. The tissue conductivities
# are physics, matched to Khadka et al., and are untouched here and asserted
# below as always. ECC is a NUMERICAL PENALTY that Workbench picks
# heuristically to enforce bonded contact: it must be large enough that the
# interface is not a resistance bottleneck, and no larger, because past that
# point it only wrecks conditioning.
#
# HOW TO CHOOSE IT -- the test is solution-independence, not convergence.
# Because ECC is numerical, the right question is not "which value solves" but
# "below which value does the answer stop changing". So this is a sweep of
# three values two orders apart, run as an array with %1 because the tlv
# account is capped at 24 CPUs and they would serialise anyway:
#
#     4e12   4 orders down, pivot ratio ~5e12
#     4e10   6 orders down, ~5e10
#     4e08   8 orders down, ~5e08
#
# Compare the reported VOLT range and the two contact node voltages across the
# three. If they agree, we are on the plateau where contact is effectively
# bonded and the lowest converging value is the right one to keep. If the
# voltages drift as ECC falls, contact resistance has started to matter and the
# value is too low -- take the highest one that still converges.
#
#SBATCH --job-name=scs_ecc
#SBATCH --array=0-2%1
#SBATCH --account=tlv
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=600G
#SBATCH --time=04:00:00
#SBATCH --output=slurm_%A_%a.out
set -uo pipefail
module load ansys/25.1
export OMP_NUM_THREADS=1
export KMP_AFFINITY=disabled

NRANK=12
PINB=${PINB:-500}          # absolute pinball radius in micrometres
ECCS=(4.e+012 4.e+010 4.e+008)
ECC=${ECCS[$SLURM_ARRAY_TASK_ID]}
SRC=$HOME/Documents/gatech/scs_model_bipolar_current.dat
WORK=${PFSDIR:-${TMPDIR:-/tmp}}/scs_ecc_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
DEST=$HOME/Documents/SCS-Modeling/ansys_training/bigdeck_ecc/out_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_ecc${ECC}
mkdir -p "$WORK" "$DEST"
echo "=== $(date) on $(hostname); PINB=-${PINB} um; ECC=${ECC}; ranks=$NRANK ==="
free -g | head -2

echo "--- patching deck (materials untouched) ---"
time sed -e 's/\r$//' \
    -e "s/^rmod,cid,6,0\./rmod,cid,6,-${PINB}./" \
    -e "s/^rmod,tid,6,0\./rmod,tid,6,-${PINB}./" \
    -e 's|^/fclean|!/fclean  ! disabled: keep the .rst|' \
    -e "s|^\*set,_maxCond,4\.e+016|*set,_maxCond,${ECC}|" \
    "$SRC" > "$WORK/scs_model.dat"

echo "--- verify patch ---"
ok=1
grep -q "^rmod,cid,6,-${PINB}\."  "$WORK/scs_model.dat" || { echo "FAIL cid pinball"; ok=0; }
grep -q "^rmod,tid,6,-${PINB}\."  "$WORK/scs_model.dat" || { echo "FAIL tid pinball"; ok=0; }
grep -q "^!/fclean"               "$WORK/scs_model.dat" || { echo "FAIL fclean"; ok=0; }
grep -q "^\*set,_maxCond,${ECC}"  "$WORK/scs_model.dat" || { echo "FAIL ECC not applied"; ok=0; }
grep -q "^\*set,_maxCond,4\.e+016" "$WORK/scs_model.dat" && { echo "FAIL some ECC still 4e16"; ok=0; }
# materials MUST still be Khadka's originals
grep -q "^MP,RSVX,1,2\.5e-13,"    "$WORK/scs_model.dat" || { echo "FAIL metal changed!"; ok=0; }
grep -q "^MP,RSVX,5,0\.05,"       "$WORK/scs_model.dat" || { echo "FAIL insulation changed!"; ok=0; }
[ "$ok" -eq 1 ] || { echo "=== aborting: patch incomplete ==="; exit 1; }
echo -n "pinball pairs patched: "; grep -c "^rmod,...,6,-${PINB}\." "$WORK/scs_model.dat"
echo -n "ECC blocks set to ${ECC}: "; grep -c "^\*set,_maxCond,${ECC}" "$WORK/scs_model.dat"
echo "materials confirmed at original Khadka values"

cat >> "$WORK/scs_model.dat" <<'POST'

/post1
set,last
nsort,volt
*get,vmx,sort,0,max
*get,vmn,sort,0,min
/com, ==== PINBALL VOLT range: min %vmn% max %vmx% ====
*get,v786,node,786,volt
*get,v3149,node,3149,volt
/com, ==== contacts: node786 %v786% V  node3149 %v3149% V ====
*get,ylead,node,786,loc,y
nsel,s,loc,y,ylead-600,ylead+600
*get,ncnt,node,0,count
/com, ==== slab nodes: %ncnt% ====
*dim,xx,array,ncnt
*dim,zz,array,ncnt
*dim,vv,array,ncnt
*vget,xx(1),node,,loc,x
*vget,zz(1),node,,loc,z
*vget,vv(1),node,,volt
*cfopen,pinball_slice,csv
*vwrite,xx(1),zz(1),vv(1)
(E16.8,',',E16.8,',',E16.8)
*cfclos
allsel,all
fini
POST

cd "$WORK"
echo "=== solve start: $(date) ==="
time ansys251 -b -dis -mpi intelmpi -np "$NRANK" -i scs_model.dat -o scs_model.out
echo "=== exit $? at $(date) ==="
echo "--- markers ---"; grep -aE "^ ==== " scs_model.out | head -20
echo "--- pivots / convergence ---"
grep -ainE "pivot|\*\*\* error|SOLUTION IS DONE|not converged|CONVERGED|NUMBER OF ITERATIONS" scs_model.out | tail -20
echo "--- results ---"; ls -la *.csv 2>/dev/null || echo "(none)"
cp -v scs_model.out pinball_slice.csv file0.err "$DEST"/ 2>/dev/null
{ echo "work: $WORK"; ls -la "$WORK" | head -30; } > "$DEST/WORK_LOCATION.txt"
ls -la "$DEST"
