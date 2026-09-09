#!/bin/bash
# Big deck attempt 2: reconditioned materials + DIRECT sparse solver forced
# OUT-OF-CORE, so it fits the 140 GB the tlv 24-CPU cap allows.
#
# What attempt 1 (job 3796753) taught us:
#   * The deck patch now verifies (CRLF stripped first -- see run_bigdeck.sbatch).
#   * PCG ran, formed all element matrices (CP 2627 s), then:
#       "The PCG solver detects that the out of balance force goes to infinity.
#        The solution has not converged. Please check for rigid body motions."
#     In an electric problem that is the signature of FLOATING REGIONS -- bodies
#     with no conduction path to the grounded component. 18,830,747 equations.
#   * MAPDL then auto-switched to EQSLV,SPARSE and the process died immediately,
#     almost certainly out of memory: the earlier successful-to-solver runs used
#     ~320 GB, and we only have 140 GB here.
#
# Why the direct solver is still worth trying:
#   The original failure was a NEGATIVE pivot (-347), which comes from the
#   2e11 metal/insulation contrast -- now reconditioned to 1e5. Floating islands
#   produce ZERO pivots, not negative ones, and MAPDL constrains those with a
#   warning rather than aborting. So direct + reconditioning may well go through
#   where PCG cannot: an iterative solver simply cannot converge on a singular
#   system, but a direct factorisation can be regularised.
#
#   BCSOPTION,,OUTOFCORE forces the sparse solver to stream its factorisation to
#   disk. /scratch has 52 TB free, so I/O space is not a concern.
#
# Fewer MPI ranks (12, not 24) deliberately: with DMP each rank duplicates some
# structures, so fewer ranks means a smaller total footprint, which is the
# binding constraint here rather than CPU time.
#
#SBATCH --job-name=scs_bigdeck2
#SBATCH --account=tlv
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=140G
#SBATCH --time=08:00:00
#SBATCH --output=slurm_%j.out
set -uo pipefail
module load ansys/25.1
export OMP_NUM_THREADS=1
export KMP_AFFINITY=disabled

NRANK=12
SRC=$HOME/Documents/gatech/scs_model_bipolar_current.dat
WORK=${PFSDIR:-${TMPDIR:-/tmp}}/scs_bigdeck2_$SLURM_JOB_ID
DEST=$HOME/Documents/SCS-Modeling/ansys_training/bigdeck_recondition/out2_$SLURM_JOB_ID
mkdir -p "$WORK" "$DEST"
echo "=== $(date) on $(hostname); work=$WORK; ranks=$NRANK ==="
df -h "$WORK" | tail -1; free -g | head -2

# NOTE the '\n' below are sed replacement escapes, NOT shell newlines. Passing a
# real newline inside the -e expression produced
# "sed: -e expression #5, char 26: unterminated `s' command".
# And do NOT put a '/' in here: the s/// delimiter is '/', so an injected
# "/config,noeldb,1" ends the replacement early and sed reports
# "unknown option to `s'". (/config,noeldb,1 is already set at the top of the
# deck anyway.)
INJECT='eqslv,sparse\nbcsoption,,outofcore'

echo "--- patching deck ---"
time sed -e 's/\r$//' \
    -e 's/^MP,RSVX,1,2\.5e-13,/MP,RSVX,1,1.0e-08,/' \
    -e 's/^MP,RSVX,5,0\.05,/MP,RSVX,5,1.0e-03,/' \
    -e 's|^/fclean|!/fclean  ! disabled: keep the .rst|' \
    -e "0,/^solve$/s//$INJECT\nsolve/" \
    "$SRC" > "$WORK/scs_model.dat"

cat >> "$WORK/scs_model.dat" <<'POST'

/post1
set,last
nsort,volt
*get,vmx,sort,0,max
*get,vmn,sort,0,min
/com, ==== BIGDECK VOLT range: min %vmn% max %vmx% ====
*get,v786,node,786,volt
*get,v3149,node,3149,volt
/com, ==== contacts: node786 %v786% V  node3149 %v3149% V ====
*get,ylead,node,786,loc,y
/com, ==== slab centred on y = %ylead% um ====
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
/show,png
/gfile,1000
/view,1,1,1,1
/auto,1
/edge,1,0
plnsol,volt
/show,close
allsel,all
fini
POST

echo "--- verify patch (abort if anything did not land) ---"
ok=1
grep -q "^MP,RSVX,1,1\.0e-08," "$WORK/scs_model.dat" || { echo "FAIL metal";      ok=0; }
grep -q "^MP,RSVX,5,1\.0e-03," "$WORK/scs_model.dat" || { echo "FAIL insulation"; ok=0; }
grep -q "^eqslv,sparse"        "$WORK/scs_model.dat" || { echo "FAIL eqslv";      ok=0; }
grep -q "^bcsoption,,outofcore" "$WORK/scs_model.dat" || { echo "FAIL bcsoption";  ok=0; }
grep -q "^!/fclean"            "$WORK/scs_model.dat" || { echo "FAIL fclean";     ok=0; }
[ "$ok" -eq 1 ] || { echo "=== aborting: deck patch incomplete ==="; exit 1; }
grep -n -A3 "^eqslv,sparse" "$WORK/scs_model.dat" | head
echo "patch verified"

cd "$WORK"
echo "=== solve start: $(date) ==="
time ansys251 -b -dis -mpi intelmpi -np "$NRANK" -i scs_model.dat -o scs_model.out
echo "=== exit $? at $(date) ==="

echo "--- markers ---";  grep -aE "^ ==== " scs_model.out | head -20
echo "--- pivot/errors/convergence ---"
grep -ainE "pivot|\*\*\* error|solution is done|not converged|out-of-core|out of core|insufficient" scs_model.out | tail -30
echo "--- results ---"; ls -la *.rst *.csv *.png 2>/dev/null || echo "(none)"

cp -v scs_model.out "$DEST"/ 2>/dev/null
cp -v bigdeck_slice.csv "$DEST"/ 2>/dev/null
for p in *.png; do [ -e "$p" ] && cp -v "$p" "$DEST"/; done
{ echo "work: $WORK"; ls -la "$WORK" | head -40; } > "$DEST/WORK_LOCATION.txt"
ls -la "$DEST"
