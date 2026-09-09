#!/bin/bash
# Truncated sub-model of the REAL anatomy, sized to fit the 140 GB that the tlv
# 24-CPU cap allows on `batch`.
#
# Why truncate. The full 20.8M-node / 18.8M-DOF model needs ~320 GB. On `batch`
# 24 CPUs buy ~140 GB (SLURM bills CPUs in proportion to memory), and even
# BCSOPTION,,OUTOFCORE was OOM killed (job 3796774). `smp` has the memory but is
# booked until ~Sep 13. So: keep the real mesh, real materials, real contacts
# and real electrodes, but solve only a box around the lead.
#
# Physics justification and its limit. A bipolar pair is a current dipole, whose
# potential falls off as ~1/r^2, so truncating at RBOX and clamping V=0 there is
# a reasonable approximation provided RBOX >> contact spacing (~4 mm here).
# RBOX = 15 mm is ~4x the spacing. The far field WILL be distorted by the
# artificial boundary -- this figure is for the near-field spread around the
# lead, NOT for absolute impedance or for anything far from the electrodes.
#
# The box grounding is also what makes the system solvable at all: the full
# model has floating (electrically isolated) bodies, which is why PCG diverged
# with "check for rigid body motions" (job 3796753). Clamping V=0 over the whole
# truncation boundary grounds every region that reaches it.
#
#SBATCH --job-name=scs_trunc
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
RBOX=${RBOX:-15000}          # half-box in micrometres (deck is /units,uMKS)
SRC=$HOME/Documents/gatech/scs_model_bipolar_current.dat
WORK=${PFSDIR:-${TMPDIR:-/tmp}}/scs_trunc_$SLURM_JOB_ID
DEST=$HOME/Documents/SCS-Modeling/ansys_training/bigdeck_truncated/out_$SLURM_JOB_ID
mkdir -p "$WORK" "$DEST"
echo "=== $(date) on $(hostname); work=$WORK; RBOX=${RBOX} um; ranks=$NRANK ==="
free -g | head -2

# Injected immediately before the deck's `solve`. Note: no '/' characters -- the
# sed s/// delimiter is '/', and an injected slash silently truncates the
# replacement ("unknown option to `s'").
read -r -d '' INJECT <<INJ || true
*get,xl,node,786,loc,x
*get,yl,node,786,loc,y
*get,zl,node,786,loc,z
rb = ${RBOX}
tolb = 800
com, ==== truncation box centre %xl% %yl% %zl% half-width %rb% um ====
allsel,all
nsel,s,loc,x,xl-rb,xl+rb
nsel,r,loc,y,yl-rb,yl+rb
nsel,r,loc,z,zl-rb,zl+rb
esln,s,1
nsle,s
cm,subn,node
cm,sube,elem
*get,nsub,node,0,count
*get,esub,elem,0,count
com, ==== submodel: %nsub% nodes %esub% elements ====
cmsel,s,subn
nsel,r,loc,x,xl-rb,xl-rb+tolb
d,all,volt,0
cmsel,s,subn
nsel,r,loc,x,xl+rb-tolb,xl+rb
d,all,volt,0
cmsel,s,subn
nsel,r,loc,y,yl-rb,yl-rb+tolb
d,all,volt,0
cmsel,s,subn
nsel,r,loc,y,yl+rb-tolb,yl+rb
d,all,volt,0
cmsel,s,subn
nsel,r,loc,z,zl-rb,zl-rb+tolb
d,all,volt,0
cmsel,s,subn
nsel,r,loc,z,zl+rb-tolb,zl+rb
d,all,volt,0
cmsel,s,sube
cmsel,s,subn
eqslv,pcg,1e-8
INJ
# collapse to a single sed replacement string with literal \n escapes
INJECT_SED=$(printf '%s' "$INJECT" | sed ':a;N;$!ba;s/\n/\\n/g')

echo "--- patching deck ---"
time sed -e 's/\r$//' \
    -e 's/^MP,RSVX,1,2\.5e-13,/MP,RSVX,1,1.0e-08,/' \
    -e 's/^MP,RSVX,5,0\.05,/MP,RSVX,5,1.0e-03,/' \
    -e 's|^/fclean|!/fclean  ! disabled: keep the .rst|' \
    -e "0,/^solve$/s//${INJECT_SED}\nsolve/" \
    "$SRC" > "$WORK/scs_model.dat"

echo "--- verify patch ---"
ok=1
grep -q "^MP,RSVX,1,1\.0e-08," "$WORK/scs_model.dat" || { echo "FAIL metal"; ok=0; }
grep -q "^MP,RSVX,5,1\.0e-03," "$WORK/scs_model.dat" || { echo "FAIL insulation"; ok=0; }
grep -q "^cm,sube,elem"        "$WORK/scs_model.dat" || { echo "FAIL truncation"; ok=0; }
grep -q "^eqslv,pcg,1e-8"      "$WORK/scs_model.dat" || { echo "FAIL solver"; ok=0; }
grep -q "^!/fclean"            "$WORK/scs_model.dat" || { echo "FAIL fclean"; ok=0; }
[ "$ok" -eq 1 ] || { echo "=== aborting: patch incomplete ==="; exit 1; }
grep -n -A4 "truncation box centre" "$WORK/scs_model.dat" | head -8
echo "patch verified"

cat >> "$WORK/scs_model.dat" <<'POST'

/post1
set,last
nsort,volt
*get,vmx,sort,0,max
*get,vmn,sort,0,min
/com, ==== TRUNC VOLT range: min %vmn% max %vmx% ====
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
*cfopen,trunc_slice,csv
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
echo "--- errors/convergence ---"
grep -ainE "pivot|\*\*\* error|solution is done|not converged|rigid body" scs_model.out | tail -20
echo "--- results ---"; ls -la *.csv *.rst 2>/dev/null || echo "(none)"
cp -v scs_model.out trunc_slice.csv "$DEST"/ 2>/dev/null
{ echo "work: $WORK"; ls -la "$WORK" | head -30; } > "$DEST/WORK_LOCATION.txt"
ls -la "$DEST"
