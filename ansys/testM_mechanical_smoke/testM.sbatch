#!/bin/bash
#SBATCH --job-name=scs_testM
#SBATCH --account=tlv
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb
#SBATCH --time=00:30:00
#SBATCH --output=slurm_%j.out
set -uo pipefail
module load ansys/25.1
cd "$SLURM_SUBMIT_DIR"
echo "=== $(date) on $(hostname) ==="
WB=/usr/local/ansys_inc/v251/aisol/.workbench
ls -la "$WB" || { echo "no .workbench launcher"; exit 1; }
# Headless Mechanical. Try plain batch first; if it aborts (the earlier
# ~/Documents/log.out crash), retry under a virtual X server if available.
echo "--- attempt 1: plain -b ---"
timeout 20m "$WB" -DSApplet -AppModeMech -b -script "$SLURM_SUBMIT_DIR/mech_smoke.py" > mech_attempt1.log 2>&1
rc=$?; echo "attempt 1 exit: $rc"; tail -40 mech_attempt1.log
if [ $rc -ne 0 ] && command -v xvfb-run >/dev/null 2>&1; then
  echo "--- attempt 2: under xvfb-run ---"
  timeout 20m xvfb-run -a "$WB" -DSApplet -AppModeMech -b -script "$SLURM_SUBMIT_DIR/mech_smoke.py" > mech_attempt2.log 2>&1
  rc=$?; echo "attempt 2 exit: $rc"; tail -40 mech_attempt2.log
fi
echo "=== files ==="; ls -la
