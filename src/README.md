# `src/` — code for the NBF_RADO-SCS study

## Current pipeline — September 30, 2026

FreeCAD geometry → Gmsh mesh → Elmer MPI potential solution → matched dorsal/
ventral field figures. Python orchestrates the compiled mesher and solver.
NEURON axon coupling is implemented as a software milestone and remains future
work for the current FEM-only poster. Remote SSH/SLURM execution of this pipeline
is planned; the completed poster runs used the local workstation.

| directory | purpose |
|---|---|
| [`freecad/`](freecad/README.md) | anatomy utilities and parametric multi-lead design |
| [`../fem/`](../fem/README.md) | Gmsh, Elmer, persistent FreeCAD visualization and Job Manager |
| [`neuron/`](neuron/README.md) | representative Aβ, Aδ and C axons and cached-field coupling |
| [`ansys/`](ansys/README.md) | historical MAPDL/SLURM experiments and shared tissue map |

Start with the [poster workflow](../docs/poster_fast_workflow.md),
[poster outline](../docs/poster_outline.md), and
[post-poster roadmap](../fem/TODO.md). Install the requirements of the component
you use; FreeCAD and Elmer are separate native installations.

## Shared tissue definitions

`ansys/tissue_map.yaml` maps every STL filename to a tissue class, an Ansys
Engineering Data material name, a conductivity, and a colour. It is the single
source of truth for "what is this body made of". `ansys/` uses it to assign
materials in Mechanical; `freecad/` reads it (as `../ansys/tissue_map.yaml`) to
colour the document so a body looks the same in FreeCAD as it does in Mechanical.
The production FEM configuration also reads this file; mesh sizing groups do
not replace material classes. Filename coverage alone does not verify spatial
coverage inside the mesh.
It lives on the Ansys side because the material values are a mirror of Mohamed's
Engineering Data and must not drift from it.

Validate it after any edit:

    python3 src/ansys/check_tissue_map.py     # every body maps to exactly one tissue

## Conventions used throughout

Measured from the geometry, not assumed — see `freecad/check_laterality.py`:

- **+X** is anatomical **left**; the midline is at x = 56.60 mm
- **+Y** is **posterior / dorsal**
- **+Z** is **rostral**
- millimetres everywhere on the FreeCAD side; the MAPDL decks are `/units,uMKS`
  (micrometres, picoamps, Tohm·µm), so watch the conversion at the boundary

RADO's `L`/`R` filename tags are unreliable — of 72 tagged bodies, 52 are
inverted. Trust geometry, not filenames.
