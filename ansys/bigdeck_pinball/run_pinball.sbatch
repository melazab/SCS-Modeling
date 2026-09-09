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
# (ansys/find_floating_bodies.py builds a proximity graph over the 245 STL
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
# Solver: PCG. With the floating regions bonded the system should no longer be
# singular, and PCG needs a fraction of the direct solver's memory, which is
# what lets this run on `batch` within the tlv 24-CPU / ~140 GB ceiling instead
# of waiting days for `smp`.
#
#SBATCH --job-name=scs_pinball
#SBATCH --account=tlv
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=140G
#SBATCH --time=06:00:00
#SBATCH --output=slurm_%j.out
set -uo pipefail
module load ansys/25.1
export OMP_NUM_THREADS=1
export KMP_AFFINITY=disabled

NRANK=12
PINB=${PINB:-500}          # absolute pinball radius in micrometres
SRC=$HOME/Documents/gatech/scs_model_bipolar_current.dat
WORK=${PFSDIR:-${TMPDIR:-/tmp}}/scs_pinball_$SLURM_JOB_ID
DEST=$HOME/Documents/SCS-Modeling/ansys_training/bigdeck_pinball/out_$SLURM_JOB_ID
mkdir -p "$WORK" "$DEST"
echo "=== $(date) on $(hostname); PINB=-${PINB} um; ranks=$NRANK ==="
free -g | head -2

echo "--- patching deck (materials untouched) ---"
time sed -e 's/\r$//' \
    -e "s/^rmod,cid,6,0\./rmod,cid,6,-${PINB}./" \
    -e "s/^rmod,tid,6,0\./rmod,tid,6,-${PINB}./" \
    -e 's|^/fclean|!/fclean  ! disabled: keep the .rst|' \
    -e "0,/^solve\$/s||eqslv,pcg,1e-8\nsolve|" \
    "$SRC" > "$WORK/scs_model.dat"

echo "--- verify patch ---"
ok=1
grep -q "^rmod,cid,6,-${PINB}\."  "$WORK/scs_model.dat" || { echo "FAIL cid pinball"; ok=0; }
grep -q "^rmod,tid,6,-${PINB}\."  "$WORK/scs_model.dat" || { echo "FAIL tid pinball"; ok=0; }
grep -q "^eqslv,pcg,1e-8"         "$WORK/scs_model.dat" || { echo "FAIL solver"; ok=0; }
grep -q "^!/fclean"               "$WORK/scs_model.dat" || { echo "FAIL fclean"; ok=0; }
# materials MUST still be Khadka's originals
grep -q "^MP,RSVX,1,2\.5e-13,"    "$WORK/scs_model.dat" || { echo "FAIL metal changed!"; ok=0; }
grep -q "^MP,RSVX,5,0\.05,"       "$WORK/scs_model.dat" || { echo "FAIL insulation changed!"; ok=0; }
[ "$ok" -eq 1 ] || { echo "=== aborting: patch incomplete ==="; exit 1; }
echo -n "pinball pairs patched: "; grep -c "^rmod,...,6,-${PINB}\." "$WORK/scs_model.dat"
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
