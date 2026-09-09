#!/bin/bash
# Truncated sub-model, attempt 3: rebuild the electrical model cleanly instead
# of trying to reconcile three overlapping coupling mechanisms.
#
# History of this line of attack:
#   3796782  spatial box worked (19.0M -> 1.88M elements) but the cut orphaned a
#            TARGE170 from its solid: "Target element ... belongs to deformable
#            target surface, while no underlying element was found."
#   3796793  deleted contacts + CPINTF: got as far as building the submodel
#            (1,701,330 nodes / 1,128,810 elements) but died on my own invalid
#            *GET,...,cp,0,maxi.
#   3796812  same, minus the bad *GET: "Specified constraint not allowed on
#            dependent degree of freedom in coupled set" x97 -- the box-face
#            D,VOLT,0 landed on CPINTF-dependent nodes.
#   3796849  CPINTF restricted to interior nodes: now CPINTF collided with the
#            deck's own 878 electrode constraint equations -- "Constraint
#            equation 2 has node 1 = 3149 ... also a dependent in coupled set".
#
# A DOF cannot be dependent in a CE and in a CP at once, and the deck ties the
# contacts with CEs while CPINTF ties tissue interfaces with CPs. So: remove BOTH
# of the deck's coupling mechanisms inside the submodel and use exactly one.
#
#   * delete CONTA174/TARGE170          (contact-based coupling)
#   * CEDELE,ALL                        (the 878 electrode equipotential CEs)
#   * FDELE,ALL,AMPS                    (the two point current injections)
#   * CPINTF,VOLT                       (single coupling mechanism, interior only)
#   * re-inject +/-1 mA DISTRIBUTED over each contact's nodes
#
# Distributing the current is what makes dropping the CEs acceptable: the metal
# is 2.5e-7 ohm*m against 0.588 ohm*m for CSF -- over 2 million times more
# conductive -- so each contact is equipotential on its own merits rather than
# by constraint. Loads (unlike constraints) are allowed on CP-dependent DOFs --
# MAPDL sums them onto the primary node.
#
# Contacts are identified by material (MAT 1 = Metal Electrode) intersected with
# a small box around each of the deck's original injection nodes, 786 and 3149.
#
#SBATCH --job-name=scs_trunc4
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
RBOX=${RBOX:-15000}      # submodel half-box, micrometres (deck is /units,uMKS)
CPTOL=${CPTOL:-300}      # CPINTF tolerance, micrometres
EBOX=${EBOX:-2000}       # half-box used to isolate each contact, micrometres
SRC=$HOME/Documents/gatech/scs_model_bipolar_current.dat
WORK=${PFSDIR:-${TMPDIR:-/tmp}}/scs_trunc4_$SLURM_JOB_ID
DEST=$HOME/Documents/SCS-Modeling/ansys_training/bigdeck_truncated/out4_$SLURM_JOB_ID
mkdir -p "$WORK" "$DEST"
echo "=== $(date) on $(hostname); RBOX=$RBOX CPTOL=$CPTOL EBOX=$EBOX ==="
free -g | head -2

read -r -d '' INJECT <<INJ || true
*get,xa,node,786,loc,x
*get,ya,node,786,loc,y
*get,za,node,786,loc,z
*get,xb,node,3149,loc,x
*get,yb,node,3149,loc,y
*get,zb,node,3149,loc,z
xm = 0.5*(xa+xb)
ym = 0.5*(ya+yb)
zm = 0.5*(za+zb)
rb = ${RBOX}
tolb = 800
cptol = ${CPTOL}
eb = ${EBOX}
/com, ==== contactA %xa% %ya% %za%  contactB %xb% %yb% %zb% ====
allsel,all
esel,s,ename,,174
esel,a,ename,,170
*get,ncont,elem,0,count
/com, ==== deleting %ncont% contact/target elements ====
edele,all
allsel,all
cedele,all
fdele,all,amps
/com, ==== deck CEs and point current loads removed ====
allsel,all
esel,s,ename,,232
esel,r,cent,x,xm-rb,xm+rb
esel,r,cent,y,ym-rb,ym+rb
esel,r,cent,z,zm-rb,zm+rb
nsle,s
cm,subn,node
cm,sube,elem
*get,nsub,node,0,count
*get,esub,elem,0,count
/com, ==== submodel %nsub% nodes %esub% elements ====
cmsel,s,subn
nsel,u,loc,x,xm-rb,xm-rb+tolb
nsel,u,loc,x,xm+rb-tolb,xm+rb
nsel,u,loc,y,ym-rb,ym-rb+tolb
nsel,u,loc,y,ym+rb-tolb,ym+rb
nsel,u,loc,z,zm-rb,zm-rb+tolb
nsel,u,loc,z,zm+rb-tolb,zm+rb
cpintf,volt,cptol
/com, ==== cpintf done ====
cmsel,s,subn
nsel,r,loc,x,xm-rb,xm-rb+tolb
d,all,volt,0
cmsel,s,subn
nsel,r,loc,x,xm+rb-tolb,xm+rb
d,all,volt,0
cmsel,s,subn
nsel,r,loc,y,ym-rb,ym-rb+tolb
d,all,volt,0
cmsel,s,subn
nsel,r,loc,y,ym+rb-tolb,ym+rb
d,all,volt,0
cmsel,s,subn
nsel,r,loc,z,zm-rb,zm-rb+tolb
d,all,volt,0
cmsel,s,subn
nsel,r,loc,z,zm+rb-tolb,zm+rb
d,all,volt,0
cmsel,s,sube
esel,r,mat,,1
cm,elece,elem
cmsel,s,elece
esel,r,cent,x,xa-eb,xa+eb
esel,r,cent,y,ya-eb,ya+eb
esel,r,cent,z,za-eb,za+eb
*get,ea,elem,0,count
nsle,s
*get,na,node,0,count
/com, ==== contact A: %ea% elements %na% nodes ====
f,all,amps,1.0e9/na
cmsel,s,elece
esel,r,cent,x,xb-eb,xb+eb
esel,r,cent,y,yb-eb,yb+eb
esel,r,cent,z,zb-eb,zb+eb
*get,eb2,elem,0,count
nsle,s
*get,nb,node,0,count
/com, ==== contact B: %eb2% elements %nb% nodes ====
f,all,amps,-1.0e9/nb
cmsel,s,sube
cmsel,s,subn
eqslv,pcg,1e-8
INJ
INJECT_SED=$(printf '%s' "$INJECT" | sed ':a;N;$!ba;s/\n/\\n/g')

echo "--- patching deck ---"
time sed -e 's/\r$//' \
    -e 's|^/fclean|!/fclean  ! disabled|' \
    -e "0,/^solve\$/s||${INJECT_SED}\nsolve|" \
    "$SRC" > "$WORK/scs_model.dat"

echo "--- verify patch ---"
ok=1
for pat in "^cedele,all" "^fdele,all,amps" \
           "^cpintf,volt,cptol" "^f,all,amps,1\.0e9" "^eqslv,pcg,1e-8"; do
  grep -q "$pat" "$WORK/scs_model.dat" || { echo "FAIL: $pat"; ok=0; }
done
[ "$ok" -eq 1 ] || { echo "=== aborting: patch incomplete ==="; exit 1; }
echo "patch verified"

cat >> "$WORK/scs_model.dat" <<'POST'

/post1
set,last
nsort,volt
*get,vmx,sort,0,max
*get,vmn,sort,0,min
/com, ==== TRUNC3 VOLT range: min %vmn% max %vmx% ====
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
echo "--- errors ---"; grep -a -A5 "\*\*\* ERROR" scs_model.out | head -25
grep -ainE "not converged|rigid body|SOLUTION IS DONE" scs_model.out | tail -8
echo "--- results ---"; ls -la *.csv 2>/dev/null || echo "(none)"
cp -v scs_model.out trunc_slice.csv file0.err "$DEST"/ 2>/dev/null
{ echo "work: $WORK"; ls -la "$WORK" | head -30; } > "$DEST/WORK_LOCATION.txt"
ls -la "$DEST"
