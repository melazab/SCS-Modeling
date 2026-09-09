#!/bin/bash
# Full 20.8M-node model on the smp partition, where there is enough memory to
# run the direct solver in core. Queued deliberately even though smp is booked
# out several days: this is the definitive full-model run, and it costs nothing
# to have it waiting.
#
# smpt10-12 are 48 cores / 2062 GB (43 GB per core), so 24 cores -- the whole
# tlv account cap -- buys ~1 TB, comfortably above the ~320 GB the earlier
# direct runs used. On `batch` the same 24 cores only buy ~140 GB, which OOM
# killed even with BCSOPTION,,OUTOFCORE (job 3796774).
#
# Changes to the deck are the same two material lines as everywhere else
# (metal 2.5e-13 -> 1.0e-08, insulation 0.05 -> 1.0e-03 Tohm*um), CR stripped
# first because the deck is CRLF, /fclean disabled, and the slice export
# appended. No solver directive: with this much memory the default distributed
# sparse solver is the right choice, and with the contrast fixed its earlier
# negative-pivot failure should not recur.
#
#SBATCH --job-name=scs_bigdeck_smp
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

SRC=$HOME/Documents/gatech/scs_model_bipolar_current.dat
WORK=${PFSDIR:-${TMPDIR:-/tmp}}/scs_bigdeck_smp_$SLURM_JOB_ID
DEST=$HOME/Documents/SCS-Modeling/ansys_training/bigdeck_recondition/outsmp_$SLURM_JOB_ID
mkdir -p "$WORK" "$DEST"
echo "=== $(date) on $(hostname); work=$WORK ==="
free -g | head -2

echo "--- patching deck ---"
time sed -e 's/\r$//' \
    -e 's|^/fclean|!/fclean  ! disabled: keep the .rst|' \
    "$SRC" > "$WORK/scs_model.dat"

cat >> "$WORK/scs_model.dat" <<'POST'

/post1
set,last
nsort,volt
*get,vmx,sort,0,max
*get,vmn,sort,0,min
/com, ==== SMP VOLT range: min %vmn% max %vmx% ====
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
*cfopen,bigdeck_slice,csv
*vwrite,xx(1),zz(1),vv(1)
(E16.8,',',E16.8,',',E16.8)
*cfclos
allsel,all
fini
POST

ok=1
grep -q "^MP,RSVX,1,1\.0e-08," "$WORK/scs_model.dat" || { echo "FAIL metal"; ok=0; }
grep -q "^MP,RSVX,5,1\.0e-03," "$WORK/scs_model.dat" || { echo "FAIL insulation"; ok=0; }
grep -q "^!/fclean"            "$WORK/scs_model.dat" || { echo "FAIL fclean"; ok=0; }
[ "$ok" -eq 1 ] || { echo "aborting: patch incomplete"; exit 1; }
echo "patch verified"

cd "$WORK"
echo "=== solve start: $(date) ==="
time ansys251 -b -dis -mpi intelmpi -np "$SLURM_NTASKS" -i scs_model.dat -o scs_model.out
echo "=== exit $? at $(date) ==="
grep -aE "^ ==== " scs_model.out | head -20
grep -ainE "pivot|\*\*\* error|solution is done" scs_model.out | tail -20
cp -v scs_model.out bigdeck_slice.csv "$DEST"/ 2>/dev/null
{ echo "work: $WORK"; ls -la "$WORK" | head -30; } > "$DEST/WORK_LOCATION.txt"
