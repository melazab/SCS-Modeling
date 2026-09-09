#!/bin/bash
# Re-run Mohamed's existing 20.8M-node Workbench deck with two changes only:
# material conditioning, and an iterative solver. Everything else -- mesh,
# boundary conditions, currents, constraint equations -- is byte-identical to
# what he built in Workbench.
#
# WHY THIS SHAPE OF JOB:
#   The tlv account has a hard group limit of cpu=24 (`sacctmgr show assoc
#   account=tlv` -> GrpTRES cpu=24,gres/gpu=1). SLURM bills CPUs in proportion
#   to requested memory, so 450G on a 2TB/384-core node implied 75 CPUs and the
#   job sat in AssocGrpCpuLimit forever. Within 24 CPUs:
#     smp    (28-42 GB/core, plenty of RAM) -> would not start until Sep 13
#     aisc*  (better ratio)                 -> not open to the tlv account
#     batch  (~6 GB/core)                   -> 24 cores = ~140 GB, starts today
#   So: batch, 24 cores, 140G. That is far less than the ~320 GB the DIRECT
#   sparse solver used -- hence the solver change below.
#
# FIX 1 - conditioning (uMKS resistivity, Tohm*um; 1 Tohm*um = 1e6 ohm*m):
#   MAT 1 metal contacts  2.5e-13 (2.5e-7 ohm*m) -> 1.0e-08 (1e-2 ohm*m)
#   MAT 5 lead insulation 0.05    (5e4   ohm*m)  -> 1.0e-03 (1e3  ohm*m)
#   Contrast 2e11 -> 1e5. The contacts are equipotential by constraint equation
#   regardless of bulk conductivity, and 1e3 ohm*m is still ~37x dura.
#
# FIX 2 - solver. The failure was:
#   "Distributed sparse solver minimum pivot = -347.2 at node 786 VOLT ...
#    An extremely large pivot ratio has been detected by the sparse solver."
#   Steady-state electric conduction is a Laplace problem: the conductance
#   matrix is symmetric positive semi-definite, so a NEGATIVE pivot is pure
#   numerical breakdown. PCG is an iterative Krylov solver -- it never pivots,
#   and it needs a small fraction of the memory of a direct factorisation.
#   Both properties are exactly what this job needs.
#
# The script tries PCG first and falls back to the direct solver in out-of-core
# mode if PCG does not converge, so an unattended overnight run still produces
# something either way.
#
#SBATCH --job-name=scs_bigdeck
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

SRC=$HOME/Documents/gatech/scs_model_bipolar_current.dat
WORK=${PFSDIR:-${TMPDIR:-/tmp}}/scs_bigdeck_$SLURM_JOB_ID
DEST=$HOME/Documents/SCS-Modeling/ansys_training/bigdeck_recondition/out_$SLURM_JOB_ID
mkdir -p "$WORK" "$DEST"
echo "=== $(date) on $(hostname); work=$WORK ==="
df -h "$WORK" | tail -1
free -g | head -2

# ---------------------------------------------------------------- patch deck
# The deck is lowercase throughout except the MP lines, and is 4.5 GB, so this
# is a single streaming pass. NEVER cat this file.
#
# GOTCHA: the deck was written by Workbench on Windows and has CRLF line
# terminators, so every line is really "solve\r". An end-anchored regex like
# /^solve$/ therefore never matches, and the injection silently does nothing --
# while the unanchored MP/fclean substitutions still work, which makes the patch
# look successful. Strip the CR first (MAPDL is happy with LF on Linux).
patch_deck () {   # $1 = solver block to inject before the first bare `solve`
  sed -e 's/\r$//' \
      -e 's/^MP,RSVX,1,2\.5e-13,/MP,RSVX,1,1.0e-08,/' \
      -e 's/^MP,RSVX,5,0\.05,/MP,RSVX,5,1.0e-03,/' \
      -e 's|^/fclean|!/fclean  ! disabled: keep the .rst|' \
      -e "0,/^solve$/s//$1\nsolve/" \
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
}

run_mapdl () {    # $1 = tag
  cd "$WORK"
  echo "=== [$1] solve start: $(date) ==="
  time ansys251 -b -dis -mpi intelmpi -np "$SLURM_NTASKS" \
       -i scs_model.dat -o "scs_model_$1.out"
  echo "=== [$1] exit $? at $(date) ==="
  echo "--- [$1] markers / errors ---"
  grep -aE "^ ==== " "scs_model_$1.out" | head -20
  grep -ainE "pivot|\*\*\* error|solution is done|not converged|failed to converge" \
       "scs_model_$1.out" | tail -20
  cp -v "scs_model_$1.out" "$DEST"/ 2>/dev/null
  [ -f bigdeck_slice.csv ] && cp -v bigdeck_slice.csv "$DEST/bigdeck_slice_$1.csv"
  for p in *.png; do [ -e "$p" ] && cp -v "$p" "$DEST/${1}_$p"; done
  ls *.rst >/dev/null 2>&1
}

# ------------------------------------------------------------- attempt 1: PCG
# eqslv,pcg,tol  -- 1e-8 is tight for a conduction problem; msave off since the
# matrix is cheap to store for SOLID232.
echo "--- patching deck for PCG ---"
time patch_deck 'eqslv,pcg,1e-8'

# Verify every edit actually landed, and abort if not. A silently-missed
# injection wastes the whole allocation running the wrong configuration.
echo "--- verify patch ---"
grep -nE "^MP,RSVX,(1|5)," "$WORK/scs_model.dat" | head
ok=1
grep -q "^MP,RSVX,1,1\.0e-08," "$WORK/scs_model.dat" || { echo "FAIL: metal not reconditioned"; ok=0; }
grep -q "^MP,RSVX,5,1\.0e-03," "$WORK/scs_model.dat" || { echo "FAIL: insulation not reconditioned"; ok=0; }
grep -q "^eqslv,pcg,1e-8"      "$WORK/scs_model.dat" || { echo "FAIL: PCG directive not injected"; ok=0; }
grep -q "^!/fclean"            "$WORK/scs_model.dat" || { echo "FAIL: /fclean not disabled"; ok=0; }
if [ "$ok" -ne 1 ]; then
  echo "=== aborting: deck patch incomplete ==="
  grep -n "^solve" "$WORK/scs_model.dat" | head -3
  exit 1
fi
echo "patch verified: materials reconditioned, PCG injected, /fclean disabled"
grep -n -B2 -A1 "^eqslv,pcg" "$WORK/scs_model.dat" | head

run_mapdl pcg
if [ -f "$WORK/bigdeck_slice.csv" ] && [ -s "$WORK/bigdeck_slice.csv" ]; then
  echo "=== PCG produced a slice export; done ==="
else
  echo "=== PCG did not produce results; falling back to direct out-of-core ==="
  rm -f "$WORK"/*.rst "$WORK"/bigdeck_slice.csv
  patch_deck 'eqslv,sparse'$'\n''bcsoption,,outofcore'
  run_mapdl direct_ooc
fi

echo "=== FINAL ==="
ls -la "$WORK"/*.rst "$WORK"/*.csv "$WORK"/*.png 2>/dev/null
{ echo "work dir: $WORK"; ls -la "$WORK" 2>/dev/null | head -40; } > "$DEST/WORK_LOCATION.txt"
ls -la "$DEST"
