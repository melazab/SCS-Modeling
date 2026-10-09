#!/bin/bash
#SBATCH -J scs-elmer-compare
#SBATCH -A tlv
#SBATCH -N 1
#SBATCH -o %x-%j.log
# Usage: sbatch -n <tasks> --mem=<size> --time=<limit> elmer_compare.sbatch <solution_bg.npz> <work_dir> [contact ...]
set -euo pipefail
SOL=$(readlink -f "$1"); WORK=$(readlink -m "$2"); shift 2
source "$HOME/scs/SCS-Modeling/fem/hpc/scs-env.sh"
NP=${SLURM_NTASKS:-1}
MEM_GB=$(( ${SLURM_MEM_PER_NODE:-8192} * 85 / 100 / 1024 ))
echo "== $(date) node $(hostname) tasks $NP budget ${MEM_GB} GB"
cd "$SCS_ROOT/fem/scripts"
ARGS=(); [[ $# -gt 0 ]] && ARGS=(--contacts "$@")
"$SCS_PY" -u resolve_stored_solution.py "$SOL" "$WORK" --cpus "$NP" --memory-gb "$MEM_GB" "${ARGS[@]}"
echo "== DONE $(date)"
