# Environment for SCS mesh/solve workers on CWRU Pioneer. Source inside a job.
# Elmer 26.2-devel 454d4a13e + Hypre 3.0.0 (same as the workstation), built by
# ~/Repos/build_elmer_pioneer.sh; Python 3.14.4 venv with the workstation's
# numpy/scipy/gmsh versions.
source "$HOME/opt/elmerfem-26.2/elmer-env.sh"
export SCS_ELMER_BIN="$HOME/opt/elmerfem-26.2/bin"
export SCS_ROOT="$HOME/scs/SCS-Modeling"
export SCS_PY="$HOME/scs/venv/bin/python"
