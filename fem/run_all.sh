#!/usr/bin/env bash
# Reproduce the whole FreeCAD-side FEM result from scratch.
#   bash fem/run_all.sh            insulating outer boundary at the canal wall
#   SCS_BACKGROUND=1 bash fem/run_all.sh   keep the bone-conductivity background
set -euo pipefail
cd "$(dirname "$0")/.."
PY=${PY:-python3}

$PY fem/scripts/verify_solver.py      # analytic sphere -- must PASS before trusting anything
$PY fem/scripts/build_mesh.py         # gmsh graded tet mesh of the stripped domain
$PY fem/scripts/assign_and_solve.py   # tissue assignment + div(sigma grad V) = 0
$PY fem/scripts/export_results.py     # VTU, tagged .msh, structured phi grid
$PY fem/scripts/analyze.py            # per-tissue |E|, comparison with Khadka 2020
$PY fem/scripts/plot_slices.py        # figure
$PY fem/scripts/sample_example.py     # worked handover to the NEURON stage
# optional independent second solver (needs apt calculix-ccx):
# $PY fem/scripts/crosscheck_ccx.py
