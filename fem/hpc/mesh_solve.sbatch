#!/bin/bash
#SBATCH -J scs-mesh-solve
#SBATCH -A tlv
#SBATCH -N 1
#SBATCH -o %x-%j.log
# Usage: sbatch -n <tasks> --mem=<size> --time=<limit> mesh_solve.sbatch <run_dir> [mesh|solve|both]
# <run_dir> holds params.json and lead/. One MPI task per core; meshing uses the
# same cores as threads. The solve budget is 85% of the SLURM allocation, never
# the node's free RAM (other jobs share the node).
set -euo pipefail
RUN=$(readlink -f "$1"); STAGE=${2:-both}
source "$HOME/scs/SCS-Modeling/fem/hpc/scs-env.sh"
NP=${SLURM_NTASKS:-1}
MEM_GB=$(( ${SLURM_MEM_PER_NODE:-8192} * 85 / 100 / 1024 ))
echo "== $(date) node $(hostname) tasks $NP solve budget ${MEM_GB} GB stage $STAGE run $RUN"
cd "$SCS_ROOT/fem/scripts"
if [[ $STAGE == mesh || $STAGE == both ]]; then
  echo "== mesh $(date)"
  OMP_NUM_THREADS=$NP "$SCS_PY" -u mesh_preview.py "$RUN" "$RUN/params.json"
fi
if [[ $STAGE == solve || $STAGE == both ]]; then
  echo "== solve $(date)"
  "$SCS_PY" -u solve_lead.py "$RUN" --lead-dir "$RUN/lead" --background --cpus "$NP" --memory-gb "$MEM_GB" \
    | tee "$RUN/solve_lead.log"
  # The same solve_report.json the Potential Visualizer writes after a local solve.
  sed -n 's/^SOLVE_JSON //p' "$RUN/solve_lead.log" > "$RUN/solve_report.json"
fi
echo "== DONE $(date)"
