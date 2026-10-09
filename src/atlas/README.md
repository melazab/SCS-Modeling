# Atlas-defined spinal tracts on the RADO cord

RADO ships one undifferentiated white-matter body, so "dorsal column" and
"spinothalamic tract" have no geometry in the model. This directory maps the
PAM50 white-matter atlas (Spinal Cord Toolbox; De Leener 2018, Lévy 2015) onto
RADO's cord. It is post-processing only: conductivities and the mesh are
untouched, and nothing under `fem/scripts/` is modified (several of those
files are hashed by the solution provenance check; they are only imported).

Outputs go to `fem/out/atlas/` (git-ignored). SCT 7.3 is installed in place at
`~/sct_7.3` without a `PATH` edit; call it as `~/sct_7.3/bin/sct_*`.

## How to run

**Interpreter:** `src/atlas/.venv`, a plain venv with numpy, scipy, nibabel and
matplotlib. `fem/.venv` will not do: it has no nibabel. One-time setup:

```bash
python3 -m venv src/atlas/.venv && src/atlas/.venv/bin/pip install numpy scipy nibabel matplotlib
```

None of these scripts needs FreeCAD or SCT on `PATH`. Only `pam50_level_match.py`
uses SCT at all, and it just reads PAM50's NIfTI files from
`~/sct_7.3/data/PAM50/template/` (the path is set in the script).
All four read RADO's original `STL_files/` through `fem/scripts/config.py` and
`inside.py`, which are imported, never modified. They do not see the T6–T9
extension.

Run in this order (each takes seconds):

| # | command | arguments | reads | writes (under `fem/out/atlas/`) |
|---|---|---|---|---|
| 1 | `src/atlas/.venv/bin/python src/atlas/voxelize_cord.py` | `--voxel` mm (default 0.25), `--margin` mm around the CSF sac (1.0), `--out` dir | white, grey, CSF STLs; aorta (left/right check) | `vox<voxel>/rado_{cord,gm,wm}_seg.nii.gz`, `rado_synth_t2.nii.gz`, `voxelize_report.json`; stops if the aorta check finds the frame mirrored |
| 2 | `src/atlas/.venv/bin/python src/atlas/make_disc_labels.py` | `--voxel` (0.25), `--naming discs\|paper` (discs), `--dir` | step 1 output, disc STLs | `vox<voxel>/rado_disc_labels.nii.gz`, `disc_labels.json` |
| 3 | `src/atlas/.venv/bin/python src/atlas/plot_voxel_check.py` | `--voxel` (0.25) | steps 1–2; lead STLs of the two poster runs in `fem/out/lead_runs/` (IDs set in the script) | `vox<voxel>/voxel_check.png` |
| 4 | `src/atlas/.venv/bin/python src/atlas/pam50_level_match.py` | `--z` RADO z of the reference section (110) | step 1 at **0.25 mm** (its centreline); RADO white/grey STLs; PAM50 cord, grey, level maps | `pam50_level_match.json`, `pam50_level_match.png` |

Status: these build and check RADO's cord images and compare RADO's
cross-section with PAM50's levels. The PAM50 → RADO **registration is not
written yet**. It waits on the level decision below and on the parametric-cord
plan in `cad-dev`.

## Frame — checked, do not re-derive

RADO STL mm are **LPS**: +x patient **left**, +y **posterior**, +z
**superior**. +x = left is established two independent ways: the descending
aorta (`Thoracic_aorta-1.STL`) is 10.8 mm toward +x of the cord, and a
right-handed frame with +y posterior and +z superior forces it. NIfTI world is
RAS = (-x, -y, z), so the voxel-to-RAS affine is a 180° rotation about z, not a
reflection. On the poster runs, lead 1 (C1–C8) is on the patient's right and
lead 2 (C9–C16) on the left, each ~2 mm from the midline.

`voxelize_cord.py` reproduces the STL volumes to 0.5 % (grey) and 0.1 % (white)
and asserts the aorta test on every run.

## Vertebral levels — names are not settled

Disc labels follow SCT's convention (value = vertebra below the disc,
C1 = 1), using RADO's disc filenames, which make the vertebrae T10–T12 (labels
17–20, ~29 mm apart). The paper says T9–T11 and the STL filenames say "T8-10";
see `src/freecad/README.md` "Vertebral levels". `--naming paper` shifts by one.

## RADO's cord is one cross-section, swept

Perpendicular sections at z = 70, 80, 110, 137 and 150 mm agree to Dice
≥ 0.987: 8.2 × 4.8 mm, 31.3 mm² cord, 9.4 mm² grey (30 %). The cord tilts
1–17° from z along its length, so z-slices are oblique; use perpendicular
sections for anything shape-related.

Against PAM50 (`pam50_level_match.py`, cords normalised to their own extents):

| PAM50 | grey fraction | grey Dice vs RADO | spinal lemniscus, one side |
|---|---:|---:|---:|
| vertebral T6–T10 (thoracic segments) | 0.17–0.21 | 0.42–0.48 | 2.2–2.5 mm² |
| vertebral T11 | 0.25 | 0.67 | 1.4 mm² |
| vertebral T12 (lumbar enlargement) | 0.37 | 0.72 | 0.8 mm² |
| RADO | 0.30 | — | (no tracts) |

RADO's grey matter is a schematic broad "H" on a flatter cord (LR/AP 1.7 vs
PAM50's 1.3–1.45); by proportion it sits nearest the lumbar enlargement.

**Consequence for registration.** Standard vertebral-level registration would
give RADO's top vertebra (ventral lead, z ≈ 137) a thoracic tract layout and its
bottom one (dorsal lead, z ≈ 80) a lumbar layout, with a spinal lemniscus about
a third the size, on a cord that is geometrically identical at both. Per-tract
dorsal-versus-ventral differences would then partly be atlas-level artefacts.
The alternative is to register a single PAM50 level onto RADO's constant
cross-section and sweep it along the centreline. Which to use is open.

## Atlas caveats

- Lévy 2015 propagated one drawn cross-section along the cord; fasciculus
  cuneatus is present in PAM50 down to T12, where it does not exist
  anatomically. Treat gracilis + cuneatus as "dorsal columns".
- The spinal lemniscus is spinothalamic **and** spinoreticular.
- RADO has no image contrast: registration fits only white/grey outlines, so the
  internal tract layout is entirely the atlas prior. Report grey-matter Dice
  next to every per-tract number.
