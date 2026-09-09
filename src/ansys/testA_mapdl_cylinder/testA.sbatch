#!/bin/bash
#SBATCH --job-name=scs_testA
#SBATCH --account=tlv
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb
#SBATCH --time=00:40:00
#SBATCH --output=slurm_%j.out
set -uo pipefail
module load ansys/25.1
cd "$SLURM_SUBMIT_DIR"
echo "=== $(date) on $(hostname) ==="

# Two failures got us here, both before MAPDL executed a single command:
#  1) `-np 4`            -> "BAD TERMINATION ... RANK 3 ... KILLED BY SIGNAL: 9"
#  2) `-np 1`            -> mpiexec bstrap_proxy failure, then
#     "srun: fatal: cpus-per-task set by two different environment variables
#      SLURM_CPUS_PER_TASK=6 != SLURM_TRES_PER_TASK=cpu=1"
# On this cluster ansys251 shells out to Intel MPI even for -np 1, and the
# SLURM 23.11+ cpus-per-task propagation conflicts with it. `-smp` forces
# shared-memory parallel with no MPI layer at all, which is the right mode for
# a model this size anyway. Unsetting the vars belts-and-braces it.
unset SLURM_CPUS_PER_TASK SLURM_TRES_PER_TASK

time ansys251 -b -smp -np "${SLURM_JOB_CPUS_PER_NODE:-4}" -i testA.dat -o testA.out
rc=$?
echo "=== ansys251 exit code: $rc ==="
echo "--- key output lines ---"
grep -E "^ ==== |\*\*\* ERROR|\*\*\* WARNING|SOLUTION IS DONE" testA.out 2>/dev/null | head -40
echo "--- png files produced ---"
ls -la *.png 2>/dev/null || echo "(no png)"
