#!/bin/bash
#SBATCH --job-name=scs_testM2
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
timeout 35m /usr/local/ansys_inc/v251/aisol/.workbench -DSApplet -AppModeMech -b \
    -script "$SLURM_SUBMIT_DIR/mech_mesh.py" > mech_m2.log 2>&1
echo "exit: $?"
cat mech_m2.log
