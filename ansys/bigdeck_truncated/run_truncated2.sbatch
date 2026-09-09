#!/bin/bash
# Truncated sub-model, attempt 2: replace contact-based coupling with direct
# VOLT coupling across interface nodes.
#
# Attempt 1 (job 3796782) successfully cut the model down -- "Out of the
# 19019112 defined elements, only 1876857 elements are selected", a 10x
# reduction that fits 140 GB -- but died with:
#
#   *** ERROR *** Target element 14603492 belongs to deformable target surface,
#   while no underlying element was found. If intended to be a deformable pair,
#   ensure that each target element has full association (shares all the nodes)
#   with the underlying solid element.
#
# i.e. the box cut through CONTA174/TARGE170 pairs and orphaned a target from
# its solid. Rather than trying to keep every pair intact across an arbitrary
# cut, delete the contact elements inside the sub-model and couple VOLT directly
# between near-coincident interface nodes with CPINTF. For a pure conduction
# problem that is equivalent to bonded electric contact, and it has a useful
# side effect: it also bonds bodies that were electrically floating, which is
# what made PCG diverge on the full model ("check for rigid body motions").
#
# Note the sed delimiter is now '|', so '/com' and other slashes are safe in the
# injected block (attempt 1 had to write "com," which MAPDL rejected with
# "COM is not a recognized SOLUTION command", losing all the diagnostics).
#
#SBATCH --job-name=scs_trunc2
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
CPTOL=${CPTOL:-200}          # CPINTF tolerance in micrometres
SRC=$HOME/Documents/gatech/scs_model_bipolar_current.dat
WORK=${PFSDIR:-${TMPDIR:-/tmp}}/scs_trunc2_$SLURM_JOB_ID
DEST=$HOME/Documents/SCS-Modeling/ansys_training/bigdeck_truncated/out2_$SLURM_JOB_ID
mkdir -p "$WORK" "$DEST"
echo "=== $(date) on $(hostname); work=$WORK; RBOX=$RBOX um; CPTOL=$CPTOL um ==="
free -g | head -2

read -r -d '' INJECT <<INJ || true
*get,xl,node,786,loc,x
*get,yl,node,786,loc,y
*get,zl,node,786,loc,z
rb = ${RBOX}
tolb = 800
cptol = ${CPTOL}
/com, ==== box centre %xl% %yl% %zl% half-width %rb% ====
allsel,all
esel,s,ename,,174
esel,a,ename,,170
*get,ncont,elem,0,count
/com, ==== deleting %ncont% contact/target elements ====
edele,all
allsel,all
esel,s,ename,,232
esel,r,cent,x,xl-rb,xl+rb
esel,r,cent,y,yl-rb,yl+rb
esel,r,cent,z,zl-rb,zl+rb
nsle,s
cm,subn,node
cm,sube,elem
*get,nsub,node,0,count
*get,esub,elem,0,count
/com, ==== submodel %nsub% nodes %esub% elements ====
cmsel,s,subn
nsel,u,loc,x,xl-rb,xl-rb+tolb
nsel,u,loc,x,xl+rb-tolb,xl+rb
nsel,u,loc,y,yl-rb,yl-rb+tolb
nsel,u,loc,y,yl+rb-tolb,yl+rb
nsel,u,loc,z,zl-rb,zl-rb+tolb
nsel,u,loc,z,zl+rb-tolb,zl+rb
cpintf,volt,cptol
/com, ==== cpintf done on interior nodes only ====
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
INJECT_SED=$(printf '%s' "$INJECT" | sed ':a;N;$!ba;s/\n/\\n/g')

echo "--- patching deck ---"
time sed -e 's/\r$//' \
    -e 's/^MP,RSVX,1,2\.5e-13,/MP,RSVX,1,1.0e-08,/' \
    -e 's/^MP,RSVX,5,0\.05,/MP,RSVX,5,1.0e-03,/' \
    -e 's|^/fclean|!/fclean  ! disabled|' \
    -e "0,/^solve\$/s||${INJECT_SED}\nsolve|" \
    "$SRC" > "$WORK/scs_model.dat"

echo "--- verify patch ---"
ok=1
grep -q "^MP,RSVX,1,1\.0e-08," "$WORK/scs_model.dat" || { echo "FAIL metal"; ok=0; }
grep -q "^MP,RSVX,5,1\.0e-03," "$WORK/scs_model.dat" || { echo "FAIL insulation"; ok=0; }
grep -q "^cpintf,volt,cptol"   "$WORK/scs_model.dat" || { echo "FAIL cpintf"; ok=0; }
grep -q "^edele,all"           "$WORK/scs_model.dat" || { echo "FAIL edele"; ok=0; }
grep -q "^eqslv,pcg,1e-8"      "$WORK/scs_model.dat" || { echo "FAIL solver"; ok=0; }
[ "$ok" -eq 1 ] || { echo "=== aborting: patch incomplete ==="; exit 1; }
grep -n -A3 "box centre" "$WORK/scs_model.dat" | head -6
echo "patch verified"

cat >> "$WORK/scs_model.dat" <<'POST'

/post1
set,last
nsort,volt
*get,vmx,sort,0,max
*get,vmn,sort,0,min
/com, ==== TRUNC2 VOLT range: min %vmn% max %vmx% ====
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
echo "--- errors ---"; grep -a -A6 "\*\*\* ERROR" scs_model.out | head -30
grep -ainE "not converged|rigid body|solution is done" scs_model.out | tail -10
echo "--- results ---"; ls -la *.csv 2>/dev/null || echo "(none)"
# always keep the per-rank error files -- the real MAPDL error lands there
cp -v scs_model.out trunc_slice.csv file0.err "$DEST"/ 2>/dev/null
{ echo "work: $WORK"; ls -la "$WORK" | head -30; } > "$DEST/WORK_LOCATION.txt"
ls -la "$DEST"
