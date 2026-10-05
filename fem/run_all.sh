#!/usr/bin/env bash
# Reproduce the whole FreeCAD-side FEM result from scratch.
#   bash fem/run_all.sh            insulating boundary of the anatomical union
#   SCS_BACKGROUND=1 bash fem/run_all.sh   keep the bone-conductivity background
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."
# Prefer fem/.venv (see fem/README.md "Setup") if it exists -- the system
# python (3.14 as of 2026-09) is too new for scipy/gmsh/pyamg wheels on this
# box and pip refuses a system-wide install (PEP 668), so requirements.txt
# needs its own venv. Falls back to plain python3 if you know what you're doing.
if [ -z "${PY:-}" ] && [ -x "$SCRIPT_DIR/.venv/bin/python" ]; then
    PY="$SCRIPT_DIR/.venv/bin/python"
fi
PY=${PY:-python3}

"$PY" fem/scripts/verify_solver.py      # analytic sphere -- must PASS before trusting anything
"$PY" fem/scripts/build_mesh.py         # gmsh graded mesh of the full anatomical domain
"$PY" fem/scripts/run_elmer.py         # tissue assignment + Elmer current-conduction solve
"$PY" fem/scripts/export_results.py     # VTU, tagged .msh, structured phi grid
"$PY" fem/scripts/analyze.py            # per-tissue |E|, comparison with Khadka 2020
"$PY" fem/scripts/plot_slices.py        # figure
"$PY" fem/scripts/sample_example.py     # worked handover to the NEURON stage
# optional independent second/third solvers:
# "$PY" fem/scripts/crosscheck_ccx.py     needs apt calculix-ccx
# "$PY" fem/scripts/crosscheck_elmer.py   needs Elmer at /opt/elmerfem (see README "Running Elmer")
