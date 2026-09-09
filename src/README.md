# `src/` — code for the NBF_RADO-SCS study

Split by which tool the code talks to, not by what stage of the pipeline it
belongs to. Anything that manipulates geometry — the STL set, the FreeCAD
document, lead construction — is under `freecad/`. Anything that talks to the
finite-element solver — MAPDL decks, SLURM job trees, result plotting — is under
`ansys/`.

| directory | runs where | what it does |
|---|---|---|
| [`freecad/`](freecad/) | your machine, under `freecadcmd` or plain `python3` | labels bodies, measures the epidural corridor, builds parametric SCS leads, colours the model |
| [`ansys/`](ansys/) | CWRU Pioneer, via `sbatch` | electric-conduction solves of the head/spine model, and post-processing |

Each directory has its own `requirements.txt` and `README.md`, because the two
halves have genuinely different dependencies: `freecad/` needs a FreeCAD install
and cannot use pip for it, while `ansys/` needs a licensed Ansys module on a
cluster and only uses Python for checking inputs and plotting outputs. Install
whichever half you actually intend to run.

## The one file both halves share

`ansys/tissue_map.yaml` maps every STL filename to a tissue class, an Ansys
Engineering Data material name, a conductivity, and a colour. It is the single
source of truth for "what is this body made of". `ansys/` uses it to assign
materials in Mechanical; `freecad/` reads it (as `../ansys/tissue_map.yaml`) to
colour the document so a body looks the same in FreeCAD as it does in Mechanical.
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
