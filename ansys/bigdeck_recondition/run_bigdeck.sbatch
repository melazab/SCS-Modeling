#!/bin/bash
# Re-run Mohamed's existing 20.8M-node Workbench deck with ONLY the material
# conditioning changed, to test whether the negative-pivot failure was caused by
# the ~2e11 conductivity contrast between the metal contacts and the lead
# insulation.
#
# The deck's physics is already correct and is left untouched:
#   d,_CM9093,volt,0.        ground on 61,742 nodes
#   f,786,amps,+1e9          +1 mA   (uMKS current unit is pA)
#   f,3149,amps,-1e9         -1 mA
#   878 x ce,next,...,volt   equipotential coupling on both contacts
#   et,*,232                 SOLID232 electric solid
#
# Changed (uMKS resistivity, Tohm*um; 1 Tohm*um = 1e6 ohm*m):
#   MAT 1 metal contacts  2.5e-13 (2.5e-7 ohm*m)  ->  1.0e-08 (1e-2 ohm*m)
#   MAT 5 lead insulation 0.05    (5e4   ohm*m)   ->  1.0e-03 (1e3  ohm*m)
# Contrast drops 2e11 -> 1e5. The contacts are held equipotential by the
# constraint equations regardless of their bulk conductivity, and 1e3 ohm*m is
# still ~37x more resistive than dura, so the physics is essentially unchanged.
#
# Also: /fclean is commented out so the .rst survives for post-processing.
#
#SBATCH --job-name=scs_bigdeck
#SBATCH --account=tlv
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=512G
#SBATCH --time=04:00:00
#SBATCH --output=slurm_%j.out
set -uo pipefail
module load ansys/25.1
export OMP_NUM_THREADS=1
export KMP_AFFINITY=disabled

SRC=$HOME/Documents/gatech/scs_model_bipolar_current.dat
WORK=${PFSDIR:-${TMPDIR:-/tmp}}/scs_bigdeck_$SLURM_JOB_ID
mkdir -p "$WORK"
echo "=== $(date) on $(hostname); work dir $WORK ==="
df -h "$WORK" | tail -1

echo "--- patching deck (streaming, original untouched) ---"
time sed -e 's/^MP,RSVX,1,2\.5e-13,/MP,RSVX,1,1.0e-08,/' \
         -e 's/^MP,RSVX,5,0\.05,/MP,RSVX,5,1.0e-03,/' \
         -e 's/^\/fclean/!\/fclean  ! disabled: keep .rst/' \
         "$SRC" > "$WORK/scs_model.dat"

echo "--- verify the patch took ---"
grep -nE "^MP,RSVX,(1|5)," "$WORK/scs_model.dat" | head
grep -cE "^!/fclean" "$WORK/scs_model.dat"

cd "$WORK"
echo "=== solve start: $(date) ==="
time ansys251 -b -dis -mpi intelmpi -np "$SLURM_NTASKS" -i scs_model.dat -o scs_model.out
rc=$?
echo "=== solve done: $(date); exit $rc ==="

echo "--- pivot / error summary ---"
grep -inE "pivot|\*\*\* error|solution is done|not converged" scs_model.out | tail -25
echo "--- results present? ---"
ls -la *.rst *.png 2>/dev/null || echo "(no rst/png)"

DEST=$HOME/Documents/SCS-Modeling/ansys_training/bigdeck_recondition/out_$SLURM_JOB_ID
mkdir -p "$DEST"
cp -v scs_model.out "$DEST"/ 2>/dev/null
cp -v *.png "$DEST"/ 2>/dev/null
# .rst can be very large; keep it on scratch and just record where it is
echo "rst location: $WORK" > "$DEST/RST_LOCATION.txt"
ls -la *.rst >> "$DEST/RST_LOCATION.txt" 2>/dev/null
echo "=== copied summary to $DEST ==="
