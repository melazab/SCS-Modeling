# `src/freecad/` — geometry side

Everything that reads or writes the STL set and the `NBF_RADO-SCS.FCStd`
document: labelling bodies, measuring the epidural space, building parametric
SCS leads, and colouring the model to match Ansys Engineering Data.

## Install

    pip install -r requirements.txt

and a FreeCAD ≥ 1.0 install providing `freecadcmd` on PATH — see the notes at the
top of `requirements.txt`, since FreeCAD itself cannot come from pip. Developed
against FreeCAD 26.3.0 / Python 3.14.4.

## Two kinds of script

| script | interpreter | what it does |
|---|---|---|
| `apply_labels.py` | `freecadcmd` | gives all 245 bodies human-readable `Label`s from `body_aliases.yaml` |
| `apply_colors.py` | `freecadcmd` | colours bodies by tissue, using `../ansys/tissue_map.yaml` |
| `measure_corridor.py` | `freecadcmd` | ray-casts the epidural mesh, writes `epidural_corridor.json` |
| `make_scs_lead.py` | `freecadcmd` | builds an *n*-contact lead swept along the measured canal centreline |
| `check_laterality.py` | `python3` | re-derives left/right from geometry and validates `body_aliases.yaml` |
| `find_floating_bodies.py` | `python3` | proximity graph over the STL surfaces; finds electrically isolated islands |

### Running the `freecadcmd` ones

`freecadcmd` owns the command line — it consumes every flag itself, and anything
after its `--pass` makes it skip the script entirely. So options go through an
environment variable instead:

    freecadcmd src/freecad/measure_corridor.py
    MAKE_LEAD_ARGS="--side ventral --z-center 110" freecadcmd src/freecad/make_scs_lead.py
    APPLY_LABELS_ARGS="--dry-run" freecadcmd src/freecad/apply_labels.py

Two more `freecadcmd` behaviours these scripts are written around: it *imports*
the script rather than running it as `__main__` (so a `__main__` guard never
fires, and the work happens at import time), and it swallows stdout (so every
script also writes a report file next to its output). It also segfaults on exit
*after* the document is saved, which is harmless.

## Why the leads have to follow a curve

The T8–T10 spine here is kyphotic, so the epidural space migrates posteriorly as
you go rostral — about 7 mm over the 43 mm length of an 8-contact lead. A lead at
fixed *y* would walk straight out of the space and through the dura. So
`measure_corridor.py` ray-casts the epidural mesh to fit a centreline *y(z)*
(cubic, 0.018 mm RMS over 98 mm) and `make_scs_lead.py` sweeps the lead along it.

Ray casting, not surface arithmetic: RADO's compartments are **hollow shells that
tile space**, each with the inner ones carved out. Measuring a channel as
`epidural.ymax − dura.ymax` is the difference of two *outer* surfaces and is not
the thickness of anything — it once produced a confident and completely wrong
"there is no dorsal epidural space" conclusion. Casting a ray along +Y and
reading consecutive entry/exit pairs gives the real numbers:

    dorsal epidural fat   2.27 – 2.35 mm
    ventral epidural fat  1.67 – 1.74 mm

so a 1.3 mm clinical lead fits on either side with no carving of the anatomy.

## Generated output

`generated_leads/<tag>/` holds one STL per contact plus one for the insulator,
already positioned in the model's own coordinates — import them without any
transform. `dorsal_z110_8c/` and `ventral_z110_8c/` are the current pair; both
verified at 100% of sampled vertices inside epidural fat and 0% in dura.
